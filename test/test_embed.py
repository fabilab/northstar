# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
from northstar import AtlasFetcher, Averages, Subsample


def test_embed_averages():
    aname = 'Darmanis_2015_nofetal'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='subsample')

    ind = [
        0, 2, 5, 8, 10, 15, 20, 25, 28, 30, 35, 38,
        40, 45, 50, 60, 70, 75, 80, 90]
    matrix = atlas[ind].copy()

    sa = Averages(
            aname,
            n_features_per_cell_type=2,
            n_features_overdispersed=5,
            n_pcs=9,
            )
    sa.fit(matrix)
    vs = sa.embed()
    assert(vs.shape[1] == 4)


def test_embed_subsample():
    aname = 'Darmanis_2015_nofetal'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='subsample')

    ind = [
        0, 2, 5, 8, 10, 15, 20, 25, 28, 30, 35, 38,
        40, 45, 50, 60, 70, 75, 80, 90]
    matrix = atlas[ind].copy()

    sa = Subsample(
            aname,
            n_features_per_cell_type=2,
            n_features_overdispersed=5,
            n_pcs=9,
            )
    sa.fit(matrix)
    vs = sa.embed()
    assert(vs.shape[1] == 4)


if __name__ == '__main__':

    test_embed_averages()
    test_embed_subsample()
