#include <iostream>
#include <vector>
#include <stdexcept>
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

    // Check for negative elements and row sums
    for (size_t i = 0; i < rows; ++i) {
        if (P[i].size() != cols) {
            std::cout << "Matrix rows have inconsistent sizes." << std::endl;
            return false;
        }
        double row_sum = 0.0;
        for (size_t j = 0; j < cols; ++j) {
            if (P[i][j] < 0) {
                std::cout << "Matrix contains negative elements." << std::endl;
                return false;
            }
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

int PeriodFinder_2n_1(const std::vector<std::vector<double>>& P, size_t i) {
    /*
    Input: stochastic matrix P of size (n, n), i - state of P (1-based index)
    Output: d - period of state i
    */
    size_t n = P.size();
    if (i > n || i < 1) {
        throw std::invalid_argument("Invalid index i");
    }
    if (!IsStochastic(P)) {
        throw std::invalid_argument("Matrix P is not stochastic");
    }

    i -= 1; // Adjust for 0-based indexing
    std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
    // Initialize M as identity matrix
    for (size_t idx = 0; idx < n; ++idx) {
        M[idx][idx] = 1.0;
    }

    int d = 0; // Initialize period

    for (size_t k = 1; k <= 2 * n; ++k) {
        if (d == 1) {
            break;
        }
        M = mat_mult(M, P);
        if (M[i][i] > 0) {
            if (d != 0) {
                d = std::gcd(d, static_cast<int>(k));
            } else {
                d = k;
            }
        }
    }

    return d;
}

int main() {
    // Define a stochastic matrix P
    std::vector<std::vector<double>> P = {
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0}
    };

    size_t state_index = 1; // Using 1-based indexing
    try {
        int period = PeriodFinder_2n_1(P, state_index);
        std::cout << "Period of state " << state_index << ": " << period << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
