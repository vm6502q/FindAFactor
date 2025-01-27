import sympy
import numpy as np

# Small semiprime example
N = 91  # 7 * 13

# Small factor base
factor_base = [2, 3, 5, 7, 11, 13]

# Example of smooth numbers found (not guaranteed squares)
smooth_numbers = [49, 77, 91]  # 7^2, 7*11, 7*13

# Factorization of each number in the factor base
factorizations = []
for num in smooth_numbers:
    factors = [sympy.factorint(num).get(p, 0) for p in factor_base]
    factorizations.append(factors)

# Convert to numpy array (exponents mod 2)
matrix = np.array(factorizations) % 2

# Perform Gaussian elimination mod 2
reduced_matrix = matrix.copy().astype(int)
rows, cols = reduced_matrix.shape
for col in range(cols):
    pivot_row = None
    for row in range(col, rows):
        if reduced_matrix[row, col] == 1:
            pivot_row = row
            break
    if pivot_row is None:
        continue
    if pivot_row != col:
        reduced_matrix[[col, pivot_row]] = reduced_matrix[[pivot_row, col]]
    for row in range(col + 1, rows):
        if reduced_matrix[row, col] == 1:
            reduced_matrix[row] = (reduced_matrix[row] + reduced_matrix[col]) % 2

# Check for a non-trivial null space solution (linear dependency)
solution_found = np.any(np.all(reduced_matrix == 0, axis=1))

# Display results
import pandas as pd

df = pd.DataFrame(matrix, columns=factor_base, index=smooth_numbers)
print(df)

solution_found
