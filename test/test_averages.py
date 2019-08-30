# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
from semiannotate import Averages


def test_constructor():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )

    sa = Averages(
            'Enge_2017',
            matrix,
            n_neighbors=2,
            threshold_neighborhood=0.8,
            n_pcs=2,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    assert(sa is not None)


def test_arguments():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )

    sa = Averages(
            'Enge_2017',
            matrix,
            n_neighbors=2,
            threshold_neighborhood=0.8,
            n_pcs=2,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    sa._check_init_arguments()
    assert(sa is not None)


def test_merge_small():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )
    n_pcs = 2
    k = 1
    threshold = 0.8

    sa = Averages(
            'Enge_2017',
            matrix,
            n_neighbors=k,
            threshold_neighborhood=threshold,
            n_pcs=n_pcs,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    sa._check_init_arguments()
    sa.fetch_atlas_if_needed()
    sa.merge_atlas_newdata()

    assert(sa is not None)


def test_neighbors_small():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )
    n_pcs = 2
    k = 1
    threshold = 0.8

    sa = Averages(
            'Enge_2017',
            matrix,
            n_neighbors=k,
            threshold_neighborhood=threshold,
            n_pcs=n_pcs,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    sa._check_init_arguments()
    sa.fetch_atlas_if_needed()
    sa.merge_atlas_newdata()
    sa.compute_neighbors()
    neis = sa.neighbors

    assert(isinstance(neis, list))
    assert(len(neis) == int(np.sum(sa.sizes)))
    for nei in neis:
        assert(isinstance(nei, list))


def test_run_small():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )
    n_pcs = 2
    k = 1
    threshold = 0.8

    sa = Averages(
            'Enge_2017',
            matrix,
            n_neighbors=k,
            threshold_neighborhood=threshold,
            n_pcs=n_pcs,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1,
            )

    sa(select_features=False)

    assert(sa.membership == ['beta', 'alpha', 'gamma'])
