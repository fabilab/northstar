# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
from semiannotate import Subsample, AtlasFetcher


def test_run_small():
    aname = 'Darmanis_2015'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='subsample')
    matrix = atlas['counts']
    cell_types = atlas['cell_types'].values
    print(cell_types)

    sa = Subsample(aname, matrix)
    sa()

    # Nobody's perfect
    assert((cell_types == sa.membership).mean() >= 0.9)

if __name__ == '__main__':

    test_run_small()

