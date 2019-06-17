# vim: fdm=indent
#author:     Fabio Zanini
#date:       17/06/19
#content:    semiknn top level module
__all__ = ['compute_graph']


import numpy as np
from sklearn.decomposition.pca import PCA


def compute_graph(matrix, sizes, n_fixed, k, threshold, metric='correlation'):
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
    '''
    neighbors = [[] for n1 in range(N - n_fixed)]

    # TODO: implement
    return neighbors
