# vim: fdm=indent
#author:     Fabio Zanini
#date:       12/08/19
#content:    Fetch cell atlas data, averaged by cell type
__all__ = ['AtlasFetcher']

import os
import numpy as np
import pandas as pd
import requests
import io
import loompy
from anndata import AnnData
import tempfile


class AtlasFetcher(object):
    '''Fetch cell atlas data from https://iosonofabio.github.io/atlas_landmarks/'''
    atlas_table = None

    def fetch_atlas_table(self):
        '''Fetch atlas table from GitHub repo'''
        url = 'https://github.com/iosonofabio/atlas_landmarks/raw/master/table.tsv'
        r = requests.get(url)
        table = pd.read_csv(io.BytesIO(r.content), sep='\t', index_col=0)

        self.atlas_table = table

    def list_atlases(self):
        '''List atlases available on GitHub repo'''
        if self.atlas_table is None:
            self.fetch_atlas_table()
        return self.atlas_table.copy()

    def fetch_atlas(self, atlas_name, kind='average', cell_types=None):
        '''Fetch an atlas from https://iosonofabio.github.io/atlas_landmarks/

        Args:
            atlas_name (str): the name of the atlas (see atlas table)
            kind (str): must be 'average' for the average expression of each
             cell type, or 'subsample' for a subsample of the atlas.
            cell_types (None or list of str): restrict to certain cell types
             within the atlas. None (default) means all cell types will be
             retained.

        Returns:
            a dictionary with two keys. 'counts' contains the gene expression
            count table. For kind == 'average', 'number_of_cells' include the
            number of cells within each cell type. For kind == 'subsample',
            'cell_types' includes the cell types of the subsample.
        '''
        if kind not in ('average', 'subsample'):
            raise ValueError('kind must be one of "average" and "subsample"')

        if self.atlas_table is None:
            self.fetch_atlas_table()

        if kind == 'subsample':
            url = self.atlas_table.at[atlas_name, 'URL_subsample']
        else:
            url = self.atlas_table.at[atlas_name, 'URL_average']

        r = requests.get(url)

        # Use a temp file, loompy has its own quirks
        fd, path = tempfile.mkstemp()
        try:
            with os.fdopen(fd, 'wb') as tmp:
                # do stuff with temp file
                tmp.write(r.content)

            with loompy.connect(path) as dsl:
                matrix = dsl.layers[''][:, :]
                features = dsl.ra['GeneName']
                cell_types_all = dsl.ca['CellType']
                if kind == 'average':
                    n_of_cells = dsl.ca['NumberOfCells']
                else:
                    cell_names = dsl.ca['CellName']

        finally:
            os.remove(path)

        # Package into anndata
        if kind == 'subsample':
            meta = pd.Series(
                data=cell_types_all,
                index=cell_names,
                )

            if cell_types is not None:
                ind_ct = meta.isin(cell_types).values.nonzero()[0]
                meta = meta.loc[ind_ct]
                matrix = matrix[:, ind_ct]

            meta.name = 'CellType'
            obs = meta.to_frame()
            obs['CellID'] = obs.index
            var = pd.Series(features, index=features)
            var.name = 'GeneName'
            var = var.to_frame()
            adata = AnnData(
                X=matrix.T,
                obs=obs,
                var=var,
                )
            adata.obs_names = adata.obs['CellID']
            adata.var_names = adata.var['GeneName']

        else:
            number_of_cells = pd.Series(
                data=n_of_cells,
                index=cell_types_all,
                )

            if cell_types is not None:
                ind_ct = number_of_cells.index.isin(cell_types).values.nonzero()[0]
                number_of_cells = number_of_cells[ind_ct]
                matrix = matrix[:, ind_ct]

            number_of_cells.name = 'NumberOfCells'
            obs = number_of_cells.to_frame()
            obs['CellID'] = obs.index
            var = pd.Series(features, index=features)
            var.name = 'GeneName'
            var = var.to_frame()
            adata = AnnData(
                X=matrix.T,
                obs=obs,
                var=var,
                )
            adata.obs_names = adata.obs['CellID']
            adata.var_names = adata.var['GeneName']

        return adata

    def fetch_multiple_atlases(
            self,
            atlas_names,
            kind='average',
            join='keep_first',
            ):
        '''Fetch and combine multiple atlases

        Args:
            atlas_names (list of str or list of dict): the names of the atlases
             (see atlas table). If a list of dict, it can be used to fetch only
             certain cell types from certain atlases. In that case, every
             element of the list must be a dict with two key-value pairs:
             'atlas_name' is the atlas name as of the atlas table, and
             'cell_types' must be a list of cell types to retain. Example:
             atlas_names=[{'atlas_name': 'Enge_2017', 'cell_tpes': ['alpha']}]
             would load the atlas Enge_2017 and only retain alpha cells.
            kind (str): must be 'average' for the average expression of each
             cell type, or 'subsample' for a subsample of the atlas.
            join (str): must be 'keep_first', 'union', or 'intersection'. This
             argument decides what to do with features that are not present
             in all atlases. 'keep_first' keeps the features in the first
             atlas and pads the other atlases with zeros, 'union' pads every
             atlas that is missing a feature and 'intersection' only keep
             features that are in all atlases.

        '''
        if kind not in ('average', 'subsample'):
            raise ValueError('kind must be one of "average" and "subsample"')

        ds = {}
        if len(atlas_names) == 0:
            return ds

        # Fetch data for all atlases
        first_atlas = None
        for atlas_name in atlas_names:
            if isinstance(atlas_name, str):
                atlas_ctypes = None
            elif isinstance(atlas_name, dict):
                atlas_ctypes = atlas_name['cell_types']
                atlas_name = atlas_name['atlas_name']
            else:
                raise ValueError('List of atlases not understood')

            if first_atlas is None:
                first_atlas = atlas_name
            ds[atlas_name] = self.fetch_atlas(
                    atlas_name,
                    kind=kind,
                    cell_types=atlas_ctypes,
                    )

        # Get features, list of all cells, etc.
        # Rename cells to ensure there are no duplicates
        cell_names = []
        cell_names_new = []
        features_list = []
        for at, d in ds.items():
            cell_names.extend(d.obs_names.values.tolist())
            cell_names_new.extend(
                    ['{:}_{:}'.format(at, x) for x in d.obs_names])
            features_list.append(d.var_names.values)
        cell_names = np.array(cell_names)
        cell_names_new = np.array(cell_names_new)

        # Deal with non-intersecting features
        if len(features_list) == 1:
            features = features_list[0]
        elif join == 'intersection':
            features_list = [set(x) for x in features_list]
            features = np.sort(set.intersection(features_list))
        elif join == 'union':
            features_list = [set(x) for x in features_list]
            features = np.sort(set.union(features_list))
        elif join == 'keep_first':
            features = np.sort(features_list[0])

        # Fill the combined dataset, including metadata
        matrix = np.empty((len(features), len(cell_names)), np.float32)
        if kind == 'average':
            n_cells_per_type = np.empty(len(cell_names), int)
        else:
            cell_types = []
        cell_dataset = []
        i = 0
        features_set = set(features)
        for at, d in ds.items():
            n = d.X.shape[0]

            # Pad missing features
            if (join == 'union') or ((join == 'keep_first') and (at != first_atlas)):
                fea_this = set(d.var_names)
                fea_missing = list(features_set - fea_this)
                fea_all = list(d.var_names) + fea_missing
                submat = np.zeros((len(fea_all), n), np.float32)
                submat[:len(fea_this)] = d.X.T

                idx = pd.Series(np.arange(len(fea_all)), index=fea_all)
                idx = idx[features].values
                submat = submat[idx]
            else:
                submat = d[:, features].X.T

            matrix[:, i: i+n] = submat

            if kind == 'average':
                n_cells_per_type[i: i+n] = d.obs['NumberOfCells'].values
            else:
                cell_types.extend(d.obs['CellType'].values.tolist())
            cell_dataset.extend([at] * d.X.shape[0])
            i += n

        cell_dataset = pd.Series(
            data=cell_dataset,
            index=cell_names_new,
            name='Dataset',
            )

        obs = cell_dataset.to_frame()
        obs['CellID'] = obs.index
        var = pd.Series(features, index=features, name='GeneName').to_frame()
        adata = AnnData(
            X=matrix.T,
            obs=obs,
            var=var,
            )
        adata.obs_names = cell_names_new
        adata.var_names = features

        if kind == 'average':
            number_of_cells = pd.Series(
                data=n_cells_per_type,
                index=cell_names_new,
                )
            adata.obs['NumberOfCells'] = number_of_cells
        else:
            cell_types = pd.Series(
                data=cell_types,
                index=cell_names_new,
                )
            adata.obs['CellType'] = cell_types

        return adata
