import numpy as np
from scipy.sparse import csr_array, diags

def div(k:int = 2, m: int = 5, dx: float = 0.2,  *args):
    n_rows = m+2
    n_cols = m+1
    
    D = csr_array((n_rows, n_cols), dtype=float)
    # D.shape
    
    neighbors = np.zeros(k)
    neighbors[0] = 1/2 - k/2
    
    for i in range(1, k):
        neighbors[i] = neighbors[i-1] + 1
    
    # vandermonde matrix over neighbors
    A = np.vander(neighbors)
    
    # first order deriv
    b = np.zeros(k)
    b[k-2] = 1
    
    # Solve the linear system
    coeffs = np.linalg.solve(A, b)
    
    i = np.arange(k//2+1, n_rows-k//2)
    j = np.arange(k)
    data = np.tile((-1)**np.arange(k), (len(i), 1))
    diags_data = np.squeeze(data.reshape(-1, 1, k))
    D = diags(diags_data, j, shape=(n_rows, n_cols))
    
    return None

def grad(k:int, m: int, dx: float,  *args):
    
    # Assertions:
    assert k >= 2, "k >= 2"
    assert k % 2 == 0, "k % 2 = 0"
    assert m >= 2*k, f"m >= {2*k} for k = {k}"
    assert len(args) < 1, "Too many input arguments."
    
    # Dimensions of G:
    if len(args) < 1:
        nodal = False
        n_rows = m + 1
        n_cols = m + 2
    else:
        nodal = True
        n_rows = m + 1
        n_cols = n_rows
    
    return None

def lap(k:int, m: int, dx: float):
    D = div(k, m, dx)
    G = grad(k, m, dx)
    
    return D * G
