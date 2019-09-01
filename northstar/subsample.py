# vim: fdm=indent
# author:     Fabio Zanini
# date:       5/08/19
# content:    Atlas subsampling
__all__ = ['Subsample']


import numpy as np
import pandas as pd
from .fetch_atlas import AtlasFetcher

try:
    from anndata import AnnData
except ImportError:
    AnnData = None


class Subsample(object):
    '''Annotate new cell types using an even subsample of an atlas'''

    def __init__(
            self,
            atlas,
            new_data,
            features=None,
            n_features_per_cell_type=30,
            n_features_overdispersed=500,
            n_pcs=20,
            n_neighbors=10,
            distance_metric='correlation',
            threshold_neighborhood=0.8,
            clustering_metric='cpm',
            resolution_parameter=0.001,
            normalize_counts=True,
            ):
        '''Prepare the model for cell annotation

        Args:
            atlas (str, list of str, or dict): cell atlas to use. If a str,
             the corresponding cell atlas from:

             https://github.com/iosonofabio/atlas_averages/blob/master/table.tsv

             is fetched (check the first column for atlas names). If a list of
             str, multiple atlases will be fetched and combined. Only features
             that are in all atlases will be kept. If you use this feature, be
             careful to not mix atlases from different species. If a dict, it
             describes a custom cell atlas and must have two entries.
             'cell_types' is a pandas Series with the cell names as
             index and the cell types to use for each cell as values.
             'counts' is a pandas.DataFrame or an anndata.AnnData structure.
             If a DataFrame, it must have features as rows and cell names as
             columns; if an AnnData, it is reversed (AnnData uses a
             different convention) and it must have the cell types as rows
             (obs_names) and the features as columns (var_names). If an AnnData,
             it will be converted into a DataFrame.

            new_data (pandas.DataFrame or anndata.AnnData): the new data to be
             clustered. If a dataframe, t must have features as rows and
             cell names as columns (as in loom files). anndata uses the opposite
             convention, so it must have cell names as rows (obs_names) and
             features as columns (var_names) and this class will transpose it.

            features (list of str or None): list of features to select after
             normalization. If None, features will be selected automatically
             based on the two next arguments, 'n_features_per_cell_type' and
             'n_features_overdispersed'. Notice that to ensure a consistent
             normalization of the atlas and the new data, feature selection
             needs to happen after normalization, so it is not recommended to
             input a pre-feature selected matrix.

            n_features_per_cell_type (int): number of features marking each
             fixed column (atlas cell type).

            n_features_overdispersed (int): number of unbiased, overdispersed
             features from the last N - n_fixed columns.

            n_pcs (int): number of principal components to keep in the PCA

            n_neighbors (int): number of neighbors in the similarity graph

            distance_metric (str): metric to use as distance. It should be a
             metric accepted by scipy.spatial.distance.cdist.

            threshold_neighborhood (float): do not consider distances larger
             than this as neighbors

            clustering_metric (str): 'cpm' (default, Cell Potts Model) or
             'modularity'. Sets the type of partition used in the clustering
             step.

            resolution_parameter (float): number between 0 and 1 that sets
             how easy it is for the clustering algorithm to make new clusters

            normalize_counts (bool): whether to renormalize the counts at the
             merging stage to make sure atlas and new data follow the same
             normalization. Be careful if you turn this off.
        '''

        self.atlas = atlas
        self.new_data = new_data
        self.features = features
        self.n_features_per_cell_type = n_features_per_cell_type
        self.n_features_overdispersed = n_features_overdispersed
        self.n_pcs = n_pcs
        self.n_neighbors = n_neighbors
        self.distance_metric = distance_metric
        self.threshold_neighborhood = threshold_neighborhood
        self.clustering_metric = clustering_metric
        self.resolution_parameter = resolution_parameter
        self.normalize_counts = normalize_counts

    def _check_init_arguments(self):
        # Custom atlas
        at = self.atlas
        if not np.isscalar(at):
            if not isinstance(self.atlas, dict):
                raise ValueError('atlas must be a dict')
            if 'counts' not in at:
                raise ValueError('atlas must have a "counts" key')
            if 'cell_types' not in at:
                raise ValueError('atlas must have a "cell_types" key')

            # The counts can be pandas.DataFrame or anndata.AnnData
            if not isinstance(at['counts'], pd.DataFrame):
                if AnnData is None:
                    raise ValueError('atlas["counts"] must be a DataFrame')
                elif not isinstance(at['counts'], AnnData):
                    raise ValueError('atlas["counts"] must be a DataFrame'
                                     ' or AnnData object')

                # AnnData uses features as columns, to transpose and convert
                at['counts'] = at['counts'].T.to_df()

            # even within AnnData, metadata colunms are pandas.DataFrame
            if not isinstance(at['cell_types'], pd.DataFrame):
                raise ValueError('atlas["cell_types"] must be a dataframe')
            if at['counts'].shape[1] != at['cell_types'].shape[0]:
                raise ValueError(
                    'atlas counts and cell_types must have the same cells')
            if (at['counts'].columns != at['cell_types'].index).any():
                raise ValueError(
                    'atlas counts and cell_types must have the same cells')

        # Make sure new data is a dataframe
        nd = self.new_data
        if not isinstance(nd, pd.DataFrame):
            if AnnData is None:
                raise ValueError('new data must be a DataFrame')
            elif not isinstance(nd, AnnData):
                raise ValueError('new_data must be a DataFrame'
                                 ' or AnnData object')

            # AnnData uses features as columns, to transpose and convert
            self.new_data = nd = nd.T.to_df()

        nf1 = self.n_features_per_cell_type
        if not isinstance(nf1, int):
            raise ValueError('n_features_per_cell_type must be an int >= 0')
        nf2 = self.n_features_overdispersed
        if not isinstance(nf1, int):
            raise ValueError('n_features_overdispersed must be an int >= 0')
        if (nf1 < 1) and (nf2 < 1):
            raise ValueError('No features selected')

    def fetch_atlas_if_needed(self):
        '''Fetch atlas(es) if needed'''

        if isinstance(self.atlas, str):
            self.atlas = AtlasFetcher().fetch_atlas(self.atlas,
                    kind='subsample')
        elif isinstance(self.atlas, list) or isinstance(self.atlas, tuple):
            self.atlas = AtlasFetcher().fetch_multiple_atlases(self.atlas,
                    kind='subsample')

    def merge_atlas_newdata(self):
        '''Merge the averaged atlas data and the new data

        This function sets the properties:
            - n_fixed: the number of cell types in the atlas
            - n_free: the number of cell types in the new data
            - cell_types: a 1D array with the cell types of the atlas
            - features_all: a 1D array with the features that were found in
            - matrix: a 2D array with the merged counts

        NOTE: is self.normalize is True, the merged count matrix is normalized
        by 1 million total counts.
        '''

        # Intersect features
        atlas_features = self.atlas['counts'].index.values
        new_data_features = self.new_data.index.values
        features = np.intersect1d(atlas_features, new_data_features)
        self.features_all = features

        # Cells, cell types, and cell numbers
        self.cell_names = self.atlas['counts'].columns.values
        self.cell_types = self.atlas['cell_types'].values
        self.n_fixed = n_fixed = self.atlas['counts'].shape[1]
        self.n_free = n_free = self.new_data.shape[1]

        # Count matrix
        L = len(features)
        N = n_fixed + n_free
        matrix = np.empty((L, N), dtype=np.float32)
        matrix[:, :n_fixed] = self.atlas['counts'].loc[features].values
        matrix[:, n_fixed:] = self.new_data.loc[features].values
        if self.normalize_counts:
            matrix *= 1e6 / matrix.sum(axis=0)
        self.matrix = matrix

    def select_features(self):
        '''Select features that define heterogeneity of the atlas and new data

        Returns:
            ndarray of feature names.
        '''
        # Shorten arg names
        matrix = self.matrix
        features = self.features
        features_all = list(self.features_all)

        if features is not None:
            ind_features = []
            for fea in features:
                ind_features.append(features_all.index(fea))
            self.features_selected = features
            self.matrix = self.matrix[ind_features]
            return

        aa = self.cell_types
        aau = list(np.unique(aa))
        n_fixed = self.n_fixed
        nf1 = self.n_features_per_cell_type
        nf2 = self.n_features_overdispersed

        ind_features = set()

        # Atlas markers
        if len(aau) > 1:
            for au in aau:
                icol1 = (aa == au).nonzero()[0]
                icol2 = (aa != au).nonzero()[0]
                ge1 = matrix[:, icol1].mean(axis=1)
                ge2 = matrix[:, icol2].mean(axis=1)
                fold_change = np.log2(ge1 + 0.1) - np.log2(ge2 + 0.1)
                markers = np.argpartition(fold_change, -nf1)[-nf1:]
                ind_features |= set(markers)

        # Unbiased on new data
        nd_mean = matrix[:, n_fixed:].mean(axis=1)
        nd_var = matrix[:, n_fixed:].var(axis=1)
        fano = (nd_var + 1e-10) / (nd_mean + 1e-10)
        overdispersed = np.argpartition(fano, -nf2)[-nf2:]
        ind_features |= set(overdispersed)

        ind_features = list(ind_features)

        self.features_selected = self.features_all[ind_features]
        self.matrix = self.matrix[ind_features]

    def compute_neighbors(self):
        '''Compute k nearest neighbors from a matrix with fixed nodes

        Returns:
            list of lists with the first k or less indices of the neighbors for
            each free column. The length is N - n_fixed. For each now, there are
            less than k entries if no other column within the distance threshold
            were found, or if N < k.


        The algorithm proceeds as follows:
        1. whiten the matrix, i.e. subtract the mean along
        the observation axis (N) and divide by the standard dev along the same axis
        2. calculate the weighted covariance matrix
        3. calculate normal PCA on that matrix
        4. calculate the distance matrix of the N1 free columns to all N columns
        5. for each free columnsm, calculate the k neighbors from the distance
        matrix, checking for threshold
        '''
        from sklearn.decomposition import PCA
        from scipy.spatial.distance import pdist, squareform

        matrix = self.matrix
        n_fixed = self.n_fixed
        k = self.n_neighbors
        n_pcs = self.n_pcs
        metric = self.distance_metric
        threshold = self.threshold_neighborhood

        # Test input arguments
        L, N = matrix.shape
        if n_fixed >= N:
            raise ValueError('n_fixed larger or equal matrix number of columns')
        if n_pcs > min(L, N):
            raise ValueError('n_pcs greater than smaller matrix dimension, those eigenvalues are zero')

        # 0. take log
        matrix = np.log10(matrix + 0.1)

        # 1. whiten
        Xnorm = ((matrix.T - matrix.mean(axis=1)) / matrix.std(axis=1, ddof=0)).T

        # take care of non-varying components
        Xnorm[np.isnan(Xnorm)] = 0

        # 2. PCA
        pca = PCA(n_components=n_pcs)
        # rvects columns are the right singular vectors
        rvects = pca.fit_transform(Xnorm.T)

        # 4. calculate distance matrix
        # rvects is N x n_pcs. Let us calculate the end that includes only the free
        # observations and then call cdist which kills the columns. The resulting
        # matrix has dimensions N1 x N
        dmat = squareform(pdist(rvects, metric=metric))
        dmax = dmat.max()

        # 5. calculate neighbors
        neighbors = [[] for n1 in range(N)]
        for i, drow in enumerate(dmat):
            nbi = neighbors[i]
            # set distance to self as a high number, to avoid self
            drow[i] = dmax + 1

            # Find largest k negative distances (k neighbors)
            ind = np.argpartition(-drow, -k)[-k:]

            # Discard the ones beyond threshold
            ind = ind[drow[ind] <= threshold]

            # Indices are not sorted within ind, we sort them to help
            # debugging even though it is not needed
            ind = ind[np.argsort(drow[ind])]

            nbi.extend(list(ind))

        self.neighbors = neighbors

    def compute_communities(self):
        '''Compute communities from a matrix with fixed nodes

        Returns:
            None, but Subsample.membership is set as an array of int with
            size N - n_fixed with the community/cluster membership of all
            columns except the first n_fixed ones.
        '''
        import inspect
        import igraph as ig
        import leidenalg

        # Check whether this version of Leiden has fixed nodes support
        opt = leidenalg.Optimiser()
        sig = inspect.getfullargspec(opt.optimise_partition)
        if 'fixed_nodes' not in sig.args:
            raise ImportError('This version of the leidenalg module does not support fixed nodes. Please update to a later (development) version')

        matrix = self.matrix
        aa = self.cell_types
        aau = list(np.unique(aa))
        n_fixed = self.n_fixed
        clustering_metric = self.clustering_metric
        resolution_parameter = self.resolution_parameter
        neighbors = self.neighbors

        L, N = matrix.shape

        # Construct graph from the lists of neighbors
        edges_d = set()
        for i, neis in enumerate(neighbors):
            for n in neis:
                edges_d.add(frozenset((i, n)))

        edges = [tuple(e) for e in edges_d]
        g = ig.Graph(n=N, edges=edges, directed=False)

        # NOTE: initial membership is singletons except for atlas nodes, which
        # get the membership they have.
        aaun = len(aau)
        initial_membership = []
        for j in range(N):
            if j < self.n_fixed:
                mb = aau.index(aa[j])
            else:
                mb = aaun + (j - n_fixed)
            initial_membership.append(mb)

        # Compute communities with semi-supervised Leiden
        if clustering_metric == 'cpm':
            partition = leidenalg.CPMVertexPartition(
                    g,
                    resolution_parameter=resolution_parameter,
                    initial_membership=initial_membership,
                    )
        elif clustering_metric == 'modularity':
            partition = leidenalg.ModularityVertexPartition(
                    g,
                    resolution_parameter=resolution_parameter,
                    initial_membership=initial_membership,
                    )
        else:
            raise ValueError(
                'clustering_metric not understood: {:}'.format(clustering_metric))

        fixed_nodes = [int(i < n_fixed) for i in range(N)]
        opt.optimise_partition(partition, fixed_nodes=fixed_nodes)
        membership = partition.membership[n_fixed:]

        # Convert the known cell types
        lstring = len(max(self.cell_types, key=len))
        self.membership = np.array(
                [str(x) for x in membership],
                dtype='U{:}'.format(lstring))
        for i, ct in enumerate(aau):
            self.membership[self.membership == str(i)] = ct

    def estimate_closest_atlas_cell_type(self):
        '''Estimate atlas cell type closest to each new cluster'''
        from scipy.spatial.distance import cdist

        matrix = self.matrix
        n_fixed = self.n_fixed
        metric = self.distance_metric
        cell_types = self.cell_types

        # Calculate averages for the new clusters
        ct_new = list(set(self.membership) - set(cell_types))
        N = len(ct_new)
        L = matrix.shape[0]
        avg_new = np.empty((L, N), np.float32)
        for i, ct in enumerate(ct_new):
            avg_new[i] = self.matrix[:, self.membership == ct].mean(axis=1)

        avg_atl = np.empty((L, len(cell_types)), np.float32)
        for i, ct in enumerate(cell_types):
            avg_atl[i] = self.matrix[:, self.membership[:n_fixed] == ct].mean(axis=1)

        # Calculate distance matrix between new and old in the high-dimensional
        # feature-selected space
        dmat = cdist(avg_new, avg_atl, metric=metric)

        # Pick the closest
        closest = np.argmin(dmat, axis=1)

        # Give it actual names
        closest = pd.Series(cell_types[closest], index=ct_new)

        return closest

    def __call__(self):
        '''Run with subsamples of the atlas

        Returns:
            None, but this instance of Subsample acquired the property
            `membership` containing the cluster memberships (cell types) of the
            columns except the first n_fixed. The first n_fixed columns are
            assumes to have distinct memberships in the range [0, n_fixed - 1].
        '''
        self._check_init_arguments()

        self.fetch_atlas_if_needed()

        self.merge_atlas_newdata()

        self.select_features()

        self.compute_neighbors()

        self.compute_communities()
