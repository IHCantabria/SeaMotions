
# Import general usage scientific libraries
from numpy import ndarray, zeros
from numpy.linalg import solve as np_solve
from scipy.special import eval_chebyt


def eval_chebyshev_2d(x: ndarray, y: ndarray, n: int, c: ndarray)->ndarray:
    sol = 0.0
    for i in range(n):
        for j in range(n):
            sol += c[i*n+j]*eval_chebyt(i, x)*eval_chebyt(j, y)

    return sol


def eval_chebyshev_2d_filter(x: ndarray, y: ndarray, c: ndarray, ncx: ndarray, ncy: ndarray)->ndarray:
    sol = 0.0
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)*eval_chebyt(ncy[i], y)

    return sol


def eval_chebyshev_3d_filter(x: ndarray, y: ndarray, z: ndarray, c: ndarray, ncx: ndarray, ncy: ndarray,
                            ncz: ndarray)->ndarray:
    sol = 0.0
    for i in range(c.shape[0]):
        sol += c[i]*eval_chebyt(ncx[i], x)*eval_chebyt(ncy[i], y)*eval_chebyt(ncz[i], z)

    return sol


def fit_chebyshev_2d(x: ndarray, y: ndarray, f: ndarray, n: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, n*n))
    for i in range(n):
        for j in range(n):
            A[:, i*n+j] = eval_chebyt(i, x)*eval_chebyt(j, y)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C


def fit_chebyshev_3d(x: ndarray, y: ndarray, z: ndarray, f: ndarray, nx: int, ny: int, nz: int)->ndarray:
    num_points = x.shape[0]
    A = zeros((num_points, nx*ny*nz))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                A[:, i*(ny*nz)+j*nz+k] = eval_chebyt(i, x)*eval_chebyt(j, y)*eval_chebyt(k, z)

    At = A.T.dot(A)
    F = A.T.dot(f)
    C = np_solve(At, F)

    return C