# vim: fdm=indent
# author:     Fabio Zanini
# date:       30/08/19
# content:    Test fetching atlases from
#             https://iosonofabio.github.io/atlas_averages/
import numpy as np
import pandas as pd
import anndata
from northstar import AtlasFetcher, subsample_atlas


def test_subsample_Darmanis_2015():
    af = AtlasFetcher()
    atlas = af.fetch_atlas('Darmanis_2015', kind='subsample')

    X = atlas['counts'].T.values
    obs = {
        'CellID': atlas['counts'].columns.values,
        'Cell Type': atlas['cell_types'].values,
        }
    var = {
        'GeneName': atlas['counts'].index.values,
        }

    atlas = anndata.AnnData(
            X=X,
            obs=obs,
            var=var,
            )

    ss = subsample_atlas(
            atlas,
            n_cells=5,
            )

    assert(isinstance(ss, anndata.AnnData))
    assert(ss.shape[0] == 40)
