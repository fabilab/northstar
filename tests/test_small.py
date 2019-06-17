# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
from semiknn import compute_graph


def test_random_small():
    N = 200
    L = 50
    matrix = np.random.rand(L, N).astype(np.float32)
    n_fixed = 3
    sizes = np.array([10, 30, 5] + [1] * (N - n_fixed))
    n_pcs = 20
    k = 5
    threshold = 0.8

    neis = compute_graph(
            matrix, sizes, n_fixed, k, threshold, n_pcs=n_pcs,
            metric='correlation',
            )

    print(neis)
