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
import tempfile


class AtlasFetcher(object):
    '''Fetch averaged cell atlas data'''
    atlas_table = None

    def fetch_atlas_table(self):
        '''Fetch atlas table from GitHub repo'''
        url = 'https://github.com/iosonofabio/atlas_averages/raw/master/table.tsv'
        r = requests.get(url)
        table = pd.read_csv(io.BytesIO(r.content), sep='\t')
        self.atlas_table = table

    def list_atlases(self):
        '''List atlases available on GitHub repo'''
        if self.atlas_table is None:
            self.fetch_atlas_table()
        return self.atlas_table.copy()

    def fetch_atlas(self, atlas_name):
        '''Fetch an atlas from GitHub repo'''
        if self.atlas_table is None:
            self.fetch_atlas_table()

        url = self.atlas_table.at[atlas_name, 'URL']
        r = requests.get(url)

        # Use a temp file, loompy has its own quirks
        fd, path = tempfile.mkstemp()
        try:
            with os.fdopen(fd, 'wb') as tmp:
                # do stuff with temp file
                tmp.write(r.content)

            with loompy.connect(path) as dsl:
                matrix = dsl.layers[''][:, :]
                cell_types = dsl.ca['CellType']
                n_of_cells = dsl.ca['NumberOfCells']
                features = dsl.ra['GeneName']

        finally:
            os.remove(path)

        return {
            'matrix': matrix,
            'cell_types': cell_types,
            'number_of_cells': n_of_cells,
            'features': features,
            }
