# sgauss_band.py
"""
Sage implementation of Modified SGAUSS Algorithm over prime field F_p.
Provides function sgauss_modified(A, b, s, L, p).
"""
from sage.all import Matrix, vector, GF, next_prime
from random import randint
from math import ceil

def random_band(n, m, L, p, nonzero=True):
    """
    Generate a random band matrix A (n x m) over GF(p) with band width L,
    ensuring nonzero entries within each row's band, and a random RHS b.

    Returns:
        - A: Sage Matrix over GF(p)
        - b: Sage vector over GF(p)
        - s: list of pivot start indices s_i for each row
    """
    field = GF(p)
    s = [randint(0, m - L) for _ in range(n)]
    A = Matrix(field, n, m)
    if nonzero:
        b = vector(field, [field.random_element(1, p) for _ in range(n)])
    else:
        b = vector(field, [field.random_element(0, p) for _ in range(n)])
        
    for i in range(n):
        start = s[i]
        for j in range(start, min(m, start + L)):
            if nonzero:
                A[i, j] = field.random_element(1, p)
            else:
                A[i, j] = field.random_element(0, p)
    return A, b, s

def sgauss_modified(A, b, s, L, p, sort_rows=False, v=False):
    """
    Perform banded Gaussian elimination over finite field F_p.

    INPUT:
        - A: n x m Matrix over F_p (Sage Matrix)
        - b: length-n vector over F_p (Sage vector)
        - s: list of start indices s_i (0-based) for each row, len(s)=n
        - L: band width
        - p: prime modulus
        - sort_rows: if True, rows will be permuted in ascending order of s before elimination

    OUTPUT:
        - x: solution vector of length m (Sage vector) if successful
        - None if failure (no pivot found or zero pivot)
    """
    n, m = A.nrows(), A.ncols()
    # Optional row sorting by s[i]
    perm = sorted(range(n), key=lambda i: s[i])
    s = [s[i] for i in perm]
    A = A.matrix_from_rows(perm)
    b = vector(GF(p), [b[i] for i in perm])
    
    
    # invperm = [0]*n
    # for idx, pi in enumerate(perm): invperm[pi] = idx    

    pivot_cols = [None] * n

    # Forward elimination
    for i in range(n):
        if v:
            print(A)
            print()
            # print(b)
            
        # Find pivot in band
        pivot = None
        for j in range(s[i], min(m, s[i] + L)):
            if A[i, j] != 0:
                pivot = j
                break
        if pivot is None:
            return None
        pivot_cols[i] = pivot

        # Normalize pivot row
        inv_val = A[i, pivot]**(-1)
        for k in range(pivot, min(m, pivot + L)):
            A[i, k] *= inv_val
        b[i] *= inv_val

        # Eliminate subsequent rows
        for r in range(i+1, n):
            if s[r] > pivot:
                break
            factor = A[r, pivot]
            if factor != 0:
                for k in range(pivot, min(m, pivot + L)):
                    A[r, k] -= factor * A[i, k]
                b[r] -= factor * b[i]
            

    # Back substitution    
    if v:
        print(A)
        # print(b)
        
    x = vector(GF(p), [0] * m)
    for i in reversed(range(n)):
        j = pivot_cols[i]
        acc = sum(A[i, k] * x[k] for k in range(j+1, min(m, j + L)))
        x[j] = b[i] - acc

    return x, pivot_cols

logp = 50
p = next_prime(2**logp)
n, L = 65536, 30
m = ceil(1.03*n)

cnt = 0
solve_fail = 0
not_correct = 0

print(n, m, L, logp)

while cnt < 100:
    if p == 2:
        A, b, s = random_band(n, m, L, p, nonzero=False)
    else:
        A, b, s = random_band(n, m, L, p, nonzero=True)
    try:
        x, pivot_cols = sgauss_modified(A, b, s, L, p)
        # if A * x != b:
            # x, pivot_cols = sgauss_modified(A, b, s, L, p, v=True)
            # not_correct += 1
    except TypeError:
        solve_fail += 1    
    cnt += 1
    print(f"\rrunning : {cnt}-th", end="", flush=True)

print()

print(solve_fail, cnt)
# if not_correct:
#     print(not_correct, "Problems ...")