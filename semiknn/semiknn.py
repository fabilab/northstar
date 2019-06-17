# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    semiknn top level module
__all__ = ['compute_graph']


import numpy as np
from sklearn.decomposition.pca import PCA


def compute_graph(matrix, sizes, n_fixed, k, threshold, n_pcs=20, metric='correlation'):
    '''Compute knn graph from a matrix with fixed nodes

    Args:
        matrix (L x N float ndarray): matrix of data. Rows are variables
        (features/genes) and columns are observations (samples/cells).

        sizes (N array of ints): number of observations (samples/cells)
        compressed in each column.

        n_fixed (int): number of columns that are fixed. These are the first
        columns of the matrix. No knn for those columns is calculated.

        k (int): number of neighbors

        threshold (float): do not consider distances larger than this as
        neighbors

        n_pcs (int): number of principal components to keep in the weighted PCA

        metric (str or function): metric to use as distance. If a string, it
        should be a metric accepted by scipy.spatial.distance.cdist. If a
        function, it should accept a (M x N)  and a (M x N1) data matrices as
        input and return a (N x N1) distance matrix. N includes both the fixed
        and the free columns, whereas N1 = N - n_fixed only includes the free
        columns.

    Returns:
        list of lists with the first k or less indices of the neighbors for
        each free column. The length is N - n_fixed. For each now, there are
        less than k entries if no other column within the distance threshold
        were found, or if N < k.


    The algorithm proceeds as follows:
    1. whiten the matrix using weights (sizes), i.e. subtract the mean along
    the observation axis (N) and divide by the standard dev along the same axis
    2. calculate the weighted covariance matrix
    3. calculate normal PCA on that matrix
    4. calculate the distance matrix of the N1 free columns to all N columns
    5. for each free columnsm, calculate the k neighbors from the distance
    matrix, checking for threshold
    '''
    from scipy.sparse.linalg import eigsh
    from scipy.spatial.distance import cdist

    # Test input arguments
    L, N = matrix.shape
    if len(sizes) != N:
        raise ValueError('Matrix and sizes dimensions do not match')
    if n_fixed >= N:
        raise ValueError('n_fixed larger or equal matrix number of columns')
    if n_pcs > min(L, N):
        raise ValueError('n_pcs greater than smaller matrix dimension, those eigenvalues are zero')

    # 1. whiten
    weights = 1.0 * sizes / sizes.sum()
    mean_w = matrix @ weights
    var_w = (matrix**2) @ weights
    std_w = np.sqrt(var_w)
    matrix_w = ((matrix.T - mean_w) / std_w).T

    # take care of non-varying components
    matrix_w[np.isnan(matrix_w)] = 0

    # 2. weighted covariance
    # This matrix has size L x L. Typically L ~ 500 << N, so the covariance
    # L x L is much smaller than N x N, hence it's fine
    cov_w = matrix_w.T @ np.diag(weights) @ matrix_w.T

    # 3. PCA
    # lvects columns are the left singular vectors L x L (gene loadings)
    evals, lvects = eigsh(cov_w, k=n_pcs)

    # calculate the right singular vectors N x N given the left singular
    # vectors and the singular values svals = np.sqrt(evals)
    # NOTE: this is true even if we truncated the PCA via n_pcs << L
    # rvects columns are the right singular vectors
    svals = np.sqrt(evals)
    rvects = matrix_w.T @ lvects @ np.diag(1.0 / svals)

    # 4. calculate distance matrix
    # rvects is N x n_pcs. Let us calculate the end that includes only the free
    # observations and then call cdist which kills the columns. The resulting
    # matrix has dimensions N1 x N
    rvects_free = rvects[n_fixed:]
    dmat = cdist(rvects_free, rvects, metric=metric)

    # 5. calculate neighbors
    neighbors = [[] for n1 in range(N - n_fixed)]
    for i, drow in enumerate(dmat):
        nbi = neighbors[i]
        # we might not need it, but it's ok fast
        if k < len(drow):
            # Find largest k negative distances (k neighbors)
            ind = np.argpartition(-drow, -k)[-k:]

            # Discard the ones beyond threshold
            ind = ind[drow[ind] <= threshold]

            # Indices are not sorted within ind, so we need to sort them
            # in descending order by distance (more efficient in the next step)
            ind = ind[np.argsort(drow[ind])][::-1]

            # Take first k, keeping multiplicities in mind
            nei_missing = k
            while (nei_missing != 0) and (len(ind)):
                j = ind[-1]
                if j >= n_fixed:
                    nbi.append(j)
                    ind = ind[:-1]
                    nei_missing -= 1
                else:
                    size = sizes[j]
                    if size >= nei_missing:
                        nbi.extend([j] * nei_missing)
                        nei_missing = 0
                    else:
                        nbi.extend([j] * size)
                        nei_missing -= size
                        ind = ind[:-1]

    return neighbors
