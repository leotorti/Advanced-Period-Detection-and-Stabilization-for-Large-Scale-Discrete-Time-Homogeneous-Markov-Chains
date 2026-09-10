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

def PeriodFinder_2n_1(P, i):
    """
    Input: stochastic matrix P of size (n, n), i - state of P (1-based index)
    Output: d - period of state i
    """
    n = P.shape[0]
    M = np.eye(n)
    d = 0  # Initialize period

    if i > n or i < 1:
        raise ValueError('Invalid index i')

    if not IsStochastic(P):
        raise ValueError('Matrix P is not stochastic')

    i -= 1  # Adjust for 0-based indexing

    for k in range(1, 2 * n):
        if d == 1:
            break
        M = np.dot(M, P)
        if M[i, i] > 0:
            d = gcd(d, k) if d != 0 else k

    return d

# Example usage:
if __name__ == "__main__":
    # Define a stochastic matrix P
    P = np.array([[0, 1, 0],
                  [0, 0, 1],
                  [1, 0, 0]])

    state_index = 1  # Using 1-based indexing
    period = PeriodFinder_2n_1(P, state_index)
    print(f"Period of state {state_index}: {period}")

