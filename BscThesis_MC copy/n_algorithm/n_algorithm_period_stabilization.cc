#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <numeric> // For std::gcd (C++17 and above)

bool IsStochastic(const std::vector<std::vector<double>>& P) {
    /*
    This function checks if the matrix P is square and stochastic.
    A matrix is stochastic if:
    1. It is square (number of rows equals number of columns)
    2. All elements are non-negative
    3. The sum of each row is 1
    */
    size_t rows = P.size();
    if (rows == 0) {
        std::cout << "Matrix is empty." << std::endl;
        return false;
    }
    size_t cols = P[0].size();
    if (rows != cols) {
        std::cout << "Matrix is not square." << std::endl;
        return false;
    }
    // Check for negative elements
    for (size_t i = 0; i < rows; ++i) {
        if (P[i].size() != cols) {
            std::cout << "Matrix rows have inconsistent sizes." << std::endl;
            return false;
        }
        for (size_t j = 0; j < cols; ++j) {
            if (P[i][j] < 0) {
                std::cout << "Matrix contains negative elements." << std::endl;
                return false;
            }
        }
    }
    // Compute row sums
    for (size_t i = 0; i < rows; ++i) {
        double row_sum = 0.0;
        for (size_t j = 0; j < cols; ++j) {
            row_sum += P[i][j];
        }
        if (std::abs(row_sum - 1.0) >= 1e-10) {
            std::cout << "Row sums are not equal to 1." << std::endl;
            return false;
        }
    }
    return true;
}

std::vector<std::vector<double>> mat_mult(const std::vector<std::vector<double>>& A,
                                          const std::vector<std::vector<double>>& B) {
    size_t n = A.size();
    std::vector<std::vector<double>> result(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t k = 0; k < n; ++k) {
            for (size_t j = 0; j < n; ++j) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

std::vector<std::vector<std::vector<double>>> PowersTransitionP(const std::vector<std::vector<double>>& P) {
    /*
    This function takes a stochastic matrix P and returns a list R such that:
    R[0] = P, R[1] = P^2, ..., R[n-1] = P^n
    */
    if (IsStochastic(P)) {
        size_t n = P.size();
        std::vector<std::vector<std::vector<double>>> R;
        R.push_back(P);  // R[0] = P
        for (size_t k = 1; k < n; ++k) {
            R.push_back(mat_mult(R[k-1], P));  // R[k] = R[k-1] * P = P^(k+1)
        }
        return R;
    } else {
        throw std::invalid_argument("Matrix P is not stochastic");
    }
}

bool i_communicates_j(const std::vector<std::vector<std::vector<double>>>& R, size_t i, size_t j, size_t n) {
    /*
    This function checks if state i communicates with state j.
    It checks if there exists k such that (P^k)[i][j] > 0 and
    if there exists h such that (P^h)[j][i] > 0.
    */
    bool flag_ij = false;
    bool flag_ji = false;
    for (size_t k = 0; k < n; ++k) {
        if (R[k][i][j] > 0 && !flag_ij) {
            flag_ij = true;
        }
        if (R[k][j][i] > 0 && !flag_ji) {
            flag_ji = true;
        }
        if (flag_ij && flag_ji) {
            break;
        }
    }
    return flag_ij && flag_ji;
}

std::vector<std::vector<int>> DividingCommunicatingClasses(const std::vector<std::vector<std::vector<double>>>& R) {
    /*
    Input: R - list of powers of P
    Output: C - matrix whose rows are the communication classes
    */
    size_t n = R[0].size();
    std::vector<std::vector<int>> C(n, std::vector<int>(n, 0));
    std::vector<bool> inserted(n, false);
    for (size_t i = 0; i < n; ++i) {
        if (std::all_of(inserted.begin(), inserted.end(), [](bool v) { return v; })) {
            break;
        }
        if (inserted[i]) {
            continue;
        }
        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                C[i][i] = 1;
                inserted[i] = true;
            }
            if (i_communicates_j(R, i, j, n)) {
                C[i][j] = 1;
                inserted[j] = true;
            }
        }
    }
    // Remove all zero rows
    std::vector<std::vector<int>> C_nonzero;
    for (size_t i = 0; i < n; ++i) {
        if (std::any_of(C[i].begin(), C[i].end(), [](int v) { return v != 0; })) {
            C_nonzero.push_back(C[i]);
        }
    }
    return C_nonzero;
}

std::pair<int, std::vector<int>> PeriodFinder_n(const std::vector<std::vector<double>>& P, size_t i) {
    /*
    Input: stochastic matrix P of size (n, n), i - state of P (1-based index)
    Output: d - period of state i, H_i - set of states communicating with i
    */
    size_t n = P.size();
    if (i > n || i < 1) {
        throw std::invalid_argument("Invalid index i");
    }
    i -= 1;  // Adjust to 0-based indexing
    std::vector<std::vector<std::vector<double>>> R = PowersTransitionP(P);
    std::vector<std::vector<int>> C = DividingCommunicatingClasses(R);
    size_t m = C.size();  // Number of communication classes
    int d = 0;  // Initialize period to be 0
    std::vector<int> H_i;
    // Find communication class of i (H_i)
    bool found = false;
    for (size_t j = 0; j < m; ++j) {
        if (C[j][i] == 1) {
            H_i = C[j];
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::runtime_error("State does not belong to any communication class.");
    }
    // Extract indices of states in H_i
    std::vector<size_t> indices;
    for (size_t idx = 0; idx < H_i.size(); ++idx) {
        if (H_i[idx] == 1) {
            indices.push_back(idx);
        }
    }
    // Compute period d
    for (size_t power = 1; power <= n; ++power) {
        if (d == 1) {
            break;
        }
        // Check if any diagonal element is non-zero in R[power - 1] restricted to indices
        bool diag_non_zero = false;
        for (size_t idx = 0; idx < indices.size(); ++idx) {
            size_t idx_i = indices[idx];
            // Diagonal element R_power[idx_i][idx_i]
            if (R[power - 1][idx_i][idx_i] != 0) {
                diag_non_zero = true;
                break;
            }
        }
        if (diag_non_zero) {
            if (d != 0) {
                d = std::gcd(d, static_cast<int>(power));
            } else {
                d = power;
            }
        }
    }
    return std::make_pair(d, H_i);
}

int main() {
    // Define a stochastic matrix P
    std::vector<std::vector<double>> P = {
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0}
    };
    size_t state_index = 1;  // Using 1-based indexing 
    try {
        auto result = PeriodFinder_n(P, state_index);
        int period = result.first;
        std::vector<int> H_i = result.second;
        std::cout << "Period of state " << state_index << ": " << period << std::endl;
        std::cout << "States communicating with state " << state_index << ": ";
        for (size_t idx = 0; idx < H_i.size(); ++idx) {
            if (H_i[idx] == 1) {
                std::cout << idx + 1 << " ";
            }
        }
        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
