# vim: fdm=indent
# author:     Fabio Zanini
# date:       30/08/19
# content:    Test fetching atlases from
#             https://iosonofabio.github.io/atlas_averages/
import numpy as np
import pandas as pd
from anndata import AnnData
from northstar import AtlasFetcher


def test_constructor():
    af = AtlasFetcher()


def test_fetch_table():
    af = AtlasFetcher()
    af.fetch_atlas_table()

    assert(af.atlas_table is not None)


def test_list_atlases():
    af = AtlasFetcher()
    table = af.list_atlases()
    assert(table is not None)


def test_fetch_Darmanis_2015():
    af = AtlasFetcher()
    atlas = af.fetch_atlas('Darmanis_2015')
    assert(isinstance(atlas, AnnData))
    assert('NumberOfCells' in atlas.obs.columns)


def test_fetch_Enge_2017():
    af = AtlasFetcher()
    atlas = af.fetch_atlas('Enge_2017')
    assert(isinstance(atlas, AnnData))
    assert('NumberOfCells' in atlas.obs.columns)

    counts = atlas[:, ['INS', 'GCG', 'PPY']]
    assert(counts['beta', 'INS'].X > counts['acinar', 'INS'].X)


def test_fetch_multiple():
    af = AtlasFetcher()
    atlas = af.fetch_multiple_atlases(['Darmanis_2015', 'Enge_2017'])
    assert(isinstance(atlas, AnnData))
    assert('NumberOfCells' in atlas.obs.columns)


def test_fetch_Darmanis_2015_subsample():
    af = AtlasFetcher()
    atlas = af.fetch_atlas('Darmanis_2015', kind='subsample')
    assert(isinstance(atlas, AnnData))
    assert('CellType' in atlas.obs.columns)


def test_fetch_multiple_subsample():
    af = AtlasFetcher()
    atlas = af.fetch_multiple_atlases(
            ['Darmanis_2015', 'Enge_2017'],
            kind='subsample')
    assert(isinstance(atlas, AnnData))
    assert('CellType' in atlas.obs.columns)
    assert('Dataset' in atlas.obs.columns)
