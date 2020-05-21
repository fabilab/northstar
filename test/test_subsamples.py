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
    matrix = atlas.copy()
    cell_types = atlas.obs['CellType'].values

    sa = Subsample(aname)
    sa.fit(matrix)

    # Small samples are a bit tricky, one cell can tip the balance
    assert((cell_types == sa.membership).mean() >= 0.85)


def test_run_across_atlas():
    atlas = AtlasFetcher().fetch_atlas(
            'Enge_2017', kind='subsample')
    matrix = atlas.copy()
    cell_types = atlas.obs['CellType'].values

    sa = Subsample(
        'Baron_2016',
        n_pcs=25,
        n_features_per_cell_type=3,
        n_features_overdispersed=200,
        )
    sa.fit(matrix)

    # Nobody's perfect
    # Baron annotates Stellate cells more accurately, so we skip them
    assert((cell_types == sa.membership)[:60].mean() >= 0.5)


if __name__ == '__main__':

    test_run_across_atlas()
