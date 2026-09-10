import numpy as np
from math import gcd

def IsStochastic(P):
    """
    This function checks if the matrix P is square and stochastic.
    A matrix is stochastic if:
    1. It is square (number of rows equals number of columns)
    2. All elements are non-negative
    3. The sum of each row is 1
    """
    rows, cols = P.shape
    if rows != cols:
        print('Matrix is not square.')
        return False

    if np.any(P < 0):
        print('Matrix contains negative elements.')
        return False

    row_sums = np.sum(P, axis=1)
    if np.all(np.abs(row_sums - 1) < 1e-10):  # Allowing for small numerical tolerance
        return True
    else:
        print('Row sums are not equal to 1.')
        return False

def PowersTransitionP(P):
    """
    This function takes a stochastic matrix P and returns a list R such that:
    R[0] = P, R[1] = P^2, ..., R[n-1] = P^n
    """
    if IsStochastic(P):
        n = P.shape[0]
        R = [P.copy()]
        for k in range(1, n):
            R.append(np.dot(R[k-1], P))  # R[k] is P^(k+1)
        return R
    else:
        raise ValueError('Matrix P is not stochastic')

def i_communicates_j(R, i, j, n):
    """
    This function checks if state i communicates with state j.
    It checks if there exists k such that (P^k)[i, j] > 0 and
    if there exists h such that (P^h)[j, i] > 0.
    """
    flag_ij = False
    flag_ji = False
    for k in range(n):
        if R[k][i, j] > 0 and not flag_ij:
            flag_ij = True
        if R[k][j, i] > 0 and not flag_ji:
            flag_ji = True
        if flag_ij and flag_ji:
            break
    return flag_ij and flag_ji

def DividingCommunicatingClasses(R):
    """
    Input: R - list of powers of P
    Output: C - matrix whose rows are the communication classes
    """
    n = R[0].shape[0]
    C = np.zeros((n, n), dtype=int)
    inserted = np.zeros(n, dtype=bool)
    for i in range(n):
        if np.all(inserted):
            break
        if inserted[i]:
            continue
        for j in range(n):
            if i == j:
                C[i, i] = 1
                inserted[i] = True
            if i_communicates_j(R, i, j, n):
                C[i, j] = 1
                inserted[j] = True
    # Remove all zero rows
    C = C[~np.all(C == 0, axis=1)]
    return C

def PeriodFinder_n(P, i):
    """
    Input: stochastic matrix P of size (n, n), i - state of P (1-based index)
    Output: d - period of state i, H_i - set of states communicating with i
    """
    n = P.shape[0]
    if i > n or i < 1:
        raise ValueError('Invalid index i')

    i -= 1  # Adjust to 0-based indexing

    R = PowersTransitionP(P)
    C = DividingCommunicatingClasses(R)
    m = C.shape[0]  # Number of communication classes
    d = 0  # Initialize period to be 0

    # Find communication class of i (H_i)
    for j in range(m):
        if C[j, i] == 1:
            H_i = C[j, :]
            break
    else:
        raise ValueError(f'State {i+1} does not belong to any communication class.')

    # Extract indices of states in H_i
    indices = np.nonzero(H_i)[0]

    # Compute period d
    for power in range(1, n + 1):
        if d == 1:
            break
        R_power = R[power - 1][np.ix_(indices, indices)]
        if np.any(np.diag(R_power) != 0):
            d = gcd(d, power) if d != 0 else power

    return d, H_i

# Example usage:
if __name__ == "__main__":
    # Define a stochastic matrix P
    P = np.array([[0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0]])

    state_index = 1  # Using 1-based indexing 
    period, H_i = PeriodFinder_n(P, state_index)
    print(f"Period of state {state_index}: {period}")
    print(f"States communicating with state {state_index}: {np.where(H_i == 1)[0] + 1}")
