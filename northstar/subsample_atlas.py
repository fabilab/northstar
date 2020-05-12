# vim: fdm=indent
# author:     Fabio Zanini
# date:       5/08/19
# content:    Atlas subsampling
__all__ = ['subsample_atlas']

import warnings
import numpy as np
import pandas as pd


def subsample_atlas(
        atlas,
        cell_type_column='CellType',
        n_cells=20,
        ):
    '''Subsample atlas across cell types

    Args:
        atlas (anndata.Anndata): The dataset to subsample.
        cell_type_column (str): The name of the metadata column to use for
            subsampling. For each unique value of this column, a certain number
            of cells will be samples (see `n_cells` argument).
        n_cells (int or dict): How many cells to sample from each cell type. If
            an int, the same number of cells will be sampled from all cell
            type. If a dict, the keys must be unique values of cell types and
            the values must be zero or positive ints that set the number of
            cells for each cell type. If a dict, the keys with zero value
            indicate what cell types to skip in the subsample.


        Note: If the atlas contains fewer than the requested number of cells
        for a cell type, all cell from that type will be sampled.
    '''
    cell_type = atlas.obs[cell_type_column]

    if np.isscalar(n_cells):
        ct_unique = np.sort(cell_type.unique())
        n_celld = {ct: n_cells for ct in ct_unique}
    else:
        n_celld = n_cells

    inds = []
    for ct, nc in n_celld.items():
        indi = (cell_type == ct).values.nonzero()[0]
        if nc < len(indi):
            indi = np.random.choice(indi, size=nc, replace=False)
            indi.sort()
        inds.extend(list(indi))

    subsample = atlas[inds]

    return subsample
