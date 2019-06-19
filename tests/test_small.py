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

    assert(isinstance(neis, list))
    assert(len(neis) == (N - n_fixed))
    for nei in neis:
        assert(isinstance(nei, list))
        assert(len(nei) <= k)


def test_adhoc_small():
    N = 200
    L = 50

    # Data matrix is random, but has defined neighborhoods
    matrix = np.ones((L, N), np.float32)
    matrix[:10, 1] += 2
    matrix[20:30, 2] -= 2
    matrix[:10, 40:80] += 2
    matrix[20:30, 100:130] -= 2

    n_fixed = 3
    sizes = np.array([10, 30, 5] + [1] * (N - n_fixed))
    n_pcs = 2
    k = 5
    threshold = 100

    neis = compute_graph(
            matrix, sizes, n_fixed, k, threshold, n_pcs=n_pcs,
            metric='euclidean',
            )

    for nei in neis[37:77]:
        assert(set(nei) & set([1] + list(range(40, 80))) == set(nei))
    for nei in neis[97:127]:
        assert(set(nei) & set([2] + list(range(100, 130))) == set(nei))

