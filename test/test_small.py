# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
from semiannotate import Averages


def test_neighbors_random():
    N = 200
    L = 50
    matrix = pd.DataFrame(index=['INS','GCG','PPY'],columns=['cell1','cell2'],data=[[2302,123],[0,5034],[0,0]])
    n_pcs = 2
    k = 1
    threshold = 0.8

    sa = Averages(
            'Darmanis_2015', matrix,
            n_neighbors=k, threshold_neighborhood=threshold, n_pcs=n_pcs,
            distance_metric='correlation',
            )
    sa.compute_neighbors()
    neis = sa.neighbors

    assert(isinstance(neis, list))
    assert(len(neis) == (N - n_fixed))
    for nei in neis:
        assert(isinstance(nei, list))
        assert(len(nei) <= k)

"""
def test_neighbors_adhoc():
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

    sa = SemiAnnotate(
            matrix, sizes=sizes, n_fixed=n_fixed,
            n_neighbors=k, threshold_neighborhood=threshold, n_pcs=n_pcs,
            distance_metric='euclidean',
            )
    sa.compute_neighbors()
    neis = sa.neighbors

    for nei in neis[37:77]:
        assert(set(nei) & set([1] + list(range(40, 80))) == set(nei))
    for nei in neis[97:127]:
        assert(set(nei) & set([2] + list(range(100, 130))) == set(nei))


def test_clustering_adhoc():
    N = 200
    L = 50

    # Data matrix is random, but has defined neighborhoods
    matrix = np.ones((L, N), np.float32)
    matrix[:10, 1] += 2
    matrix[20:30, 2] -= 2
    matrix[:10, 40:45] += 2
    matrix[20:30, 100:105] -= 2

    n_fixed = 3
    sizes = np.array([10, 30, 5] + [1] * (N - n_fixed))
    n_pcs = 2
    k = 8
    threshold = 100

    sa = SemiAnnotate(
            matrix, sizes=sizes, n_fixed=n_fixed,
            n_neighbors=k, threshold_neighborhood=threshold, n_pcs=n_pcs,
            distance_metric='euclidean',
            clustering_metric='cpm',
            resolution_parameter=0.0001,
            )
    sa.compute_neighbors()
    sa.compute_communities()
    communities = sa.membership

    assert(len(communities) == N - n_fixed)
    assert(len(set(communities[37:42])) == 1)
    assert(len(set(communities[97:102])) == 1)
    assert(communities[37] == 1)
    assert(communities[97] == 2)
    """
