# vim: fdm=indent
# author:     Fabio Zanini
# date:       5/08/19
# content:    Atlas subsampling
__all__ = ['subsample_atlas', 'average_atlas']

import warnings
import numpy as np
import pandas as pd
import anndata


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
        n_celld = cell_type.value_counts()
        ct_unique = n_celld.index.values
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


def average_atlas(
        atlas,
        cell_type_column='CellType',
        ):
    '''Subsample atlas across cell types

    Args:
        atlas (anndata.Anndata): The dataset to subsample.
        cell_type_column (str): The name of the metadata column to use for
            subsampling. For each unique value of this column, a certain number
            of cells will be samples (see `n_cells` argument).

        Note: If the atlas contains fewer than the requested number of cells
        for a cell type, all cell from that type will be sampled.
    '''
    cell_type = atlas.obs[cell_type_column]
    n_celld = cell_type.value_counts()
    ct_unique = list(n_celld.index)

    matrix = np.zeros((len(ct_unique), atlas.shape[1]), np.float32)

    for ict, ct in enumerate(ct_unique):
        indi = (cell_type == ct).values.nonzero()[0]
        mati = atlas.X[indi]
        if not isinstance(mati, np.ndarray):
            # Make a dense matrix, it's ok since there are not many types
            mati = mati.toarray()
        matrix[ict] = mati.mean(axis=0)

    ave = anndata.AnnData(
            X=matrix,
            obs={'CellType': ct_unique, 'NumberOfCells': n_celld.values},
            )
    ave.obs_names = ct_unique
    ave.var = atlas.var

    return ave
