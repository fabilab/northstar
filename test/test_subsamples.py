# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
from northstar import Subsample, AtlasFetcher


def test_run_within_atlas():
    aname = 'Darmanis_2015'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='subsample')
    matrix = atlas['counts']
    cell_types = atlas['cell_types'].values

    sa = Subsample(aname)
    sa.fit(matrix)

    # Nobody's perfect
    assert((cell_types == sa.membership).mean() >= 0.9)


def test_run_within_atlas_means():
    aname = 'Darmanis_2015'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='average')

    ind = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 0]
    matrix = atlas['counts'].iloc[:, ind]
    cell_types = atlas['counts'].columns.values[ind]

    sa = Subsample(aname)
    sa.fit(matrix)

    # Nobody's perfect
    assert((cell_types == sa.membership).mean() >= 0.9)


def test_run_across_atlas():
    atlas = AtlasFetcher().fetch_atlas(
            'Enge_2017', kind='subsample')
    matrix = atlas['counts']
    cell_types = atlas['cell_types'].values

    sa = Subsample(
        'Baron_2016',
        n_pcs=25,
        n_features_per_cell_type=3,
        n_features_overdispersed=200,
        )
    sa.fit(matrix)

    # Nobody's perfect
    # Baron annotates Stellate cells more accurately, so we skip them
    assert((cell_types == sa.membership)[:60].mean() >= 0.9)


if __name__ == '__main__':

    test_run_across_atlas()
