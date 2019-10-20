# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Atlas averages
__all__ = ['Averages']


import numpy as np
import pandas as pd
from .fetch_atlas import AtlasFetcher

try:
    from anndata import AnnData
except ImportError:
    AnnData = None


class Averages(object):
    '''Annotate new cell types using averages of an atlas'''

    def __init__(
            self,
            atlas,
            n_cells_per_type=None,
            features=None,
            n_features_per_cell_type=30,
            n_features_overdispersed=500,
            features_additional=None,
            n_pcs=20,
            n_neighbors=10,
            n_neighbors_out_of_atlas=5,
            distance_metric='correlation',
            threshold_neighborhood=0.8,
            clustering_metric='cpm',
            resolution_parameter=0.001,
            normalize_counts=True,
            join='keep_first',
            ):
        '''Prepare the model for cell annotation

        Args:
            atlas (str, list of str, list of dict, or dict): cell atlas to use.
            If a str, the corresponding cell atlas from:

             https://github.com/iosonofabio/atlas_averages/blob/master/table.tsv

             is fetched (check the first column for atlas names). If a list of
             str, multiple atlases will be fetched and combined. Only features
             that are in all atlases will be kept. If you use this feature, be
             careful to not mix atlases from different species. If a list of
             dict, it merges atlases as above but you can specify what cell
             types to fetch from each atlas. Each element of the list must be a
             dict with two key-value pairs: 'atlas_name' is the atlas name, and
             'cell_types' must be a list of cell types to retain. Example:
             atlas=[{'atlas_name': 'Enge_2017', 'cell_tpes': ['alpha']}] would
             load the atlas Enge_2017 and only retain alpha cells. If a dict,
             it can refer to two options. The first is a single atlas with a
             specification to retain only certain cell types. The format is as
             above, e.g. to select only alpha cells from Enge_2017 you can use:
             atlas={'atlas_name': 'Enge_2017', 'cell_tpes': ['alpha']}.
             The second option describes a custom cell atlas. In this case, the
             dict must have two entries, 'number_of_cells' and 'counts'.
             'number_of_cells' is a pandas Series with the cell types as index
             and the number of cells to use for each cell type as values.
             'counts' is a pandas.DataFrame or an anndata.AnnData structure.
             If a DataFrame, it must have features as rows and cell types as
             columns; if an AnnData, it is reversed (AnnData uses a
             different convention) and it must have the cell types as rows
             (obs_names) and the features as columns (var_names). If an AnnData,
             it will be converted into a DataFrame.

            n_cells_per_type (None or int): if None, use the number of cells
             per type from the atlas. Else, fix it to this number for all types.

            features (list of str or None): list of features to select after
             normalization. If None, features will be selected automatically
             based on the two next arguments, 'n_features_per_cell_type' and
             'n_features_overdispersed'. Notice that to ensure a consistent
             normalization of the atlas and the new data, feature selection
             needs to happen after normalization, so it is not recommended to
             input a pre-feature selected matrix.

            n_features_per_cell_type (int): number of features marking each
             fixed column (atlas cell type). The argument 'features' takes
             priority over this one.

            n_features_overdispersed (int): number of unbiased, overdispersed
             features to be picked from the new dataset. The argument
             'features' takes priority over this one.

            features_additional (list of str or None): additional features to
             keep on top of automatic selection. The argument 'features' takes
             priority over this one.

            n_pcs (int): number of principal components to keep in the weighted
             PCA.

            n_neighbors (int): number of neighbors in the similarity graph.

            n_neighbors_out_of_atlas (int): number of neighbors coming out of
             the atlas nodes into the new dataset.

            distance_metric (str): metric to use as distance. It should be a
             metric accepted by scipy.spatial.distance.cdist.

            threshold_neighborhood (float): do not consider distances larger than this as
             neighbors

            clustering_metric (str): 'cpm' (default, Cell Potts Model) or
             'modularity'. Sets the type of partition used in the clustering
             step.

            resolution_parameter (float): number between 0 and 1 that sets
             how easy it is for the clustering algorithm to make new clusters

            normalize_counts (bool): whether to renormalize the counts at the
             merging stage to make sure atlas and new data follow the same
             normalization. Be careful if you turn this off.

            join (str): must be 'keep_first', 'union', or 'intersection'. This
             argument is used when sourcing multiple atlases and decides what
             to do with features that are not present in all atlases.
             'keep_first' keeps the features in the first atlas and pads the
             other atlases with zeros, 'union' pads every atlas that is missing
             a feature and 'intersection' only keep features that are in all
             atlases.
        '''

        self.atlas = atlas
        self.n_cells_per_type = n_cells_per_type
        self.features = features
        self.n_features_per_cell_type = n_features_per_cell_type
        self.n_features_overdispersed = n_features_overdispersed
        self.features_additional = features_additional
        self.n_pcs = n_pcs
        self.n_neighbors = n_neighbors
        self.n_neighbors_out_of_atlas = n_neighbors_out_of_atlas
        self.distance_metric = distance_metric
        self.threshold_neighborhood = threshold_neighborhood
        self.clustering_metric = clustering_metric
        self.resolution_parameter = resolution_parameter
        self.normalize_counts = normalize_counts
        self.join = join

    def fit(self, new_data):
        '''Run with averages of the atlas

        Args:
            new_data (pandas.DataFrame or anndata.AnnData): the new data to be
             clustered. If a dataframe, t must have features as rows and
             cell names as columns (as in loom files). anndata uses the opposite
             convention, so it must have cell names as rows (obs_names) and
             features as columns (var_names) and this class will transpose it.

        Returns:
            None, but this instance of Averages acquired the property
            `membership` containing the cluster memberships (cell types) of the
            columns except the first n_fixed. The first n_fixed columns are
            assumes to have distinct memberships in the range [0, n_fixed - 1].
        '''
        self.new_data = new_data

        self._check_init_arguments()
        self.fetch_atlas_if_needed()
        self.merge_atlas_newdata()
        self.select_features()
        self.compute_pca()
        self.compute_similarity_graph()
        self.cluster_graph()

    def fit_transform(self, new_data):
        '''Run with averages of the atlas and return the cell types

        Args:
            new_data (pandas.DataFrame or anndata.AnnData): the new data to be
             clustered. If a dataframe, t must have features as rows and
             cell names as columns (as in loom files). anndata uses the opposite
             convention, so it must have cell names as rows (obs_names) and
             features as columns (var_names) and this class will transpose it.

        Returns:
            the cluster memberships (cell types) of the
            columns except the first n_fixed. The first n_fixed columns are
            assumes to have distinct memberships in the range [0, n_fixed - 1].
        '''
        self.fit(new_data)
        return self.membership

    def _check_init_arguments(self):
        # Custom atlas
        at = self.atlas
        if isinstance(at, str):
            pass
        elif isinstance(at, list) or isinstance(at, tuple):
            for elem in at:
                if isinstance(elem, str):
                    pass
                elif isinstance(elem, dict):
                    if 'atlas_name' not in elem:
                        raise ValueError('List of atlases: format incorrect')
                    if 'cell_types' not in elem:
                        raise ValueError('List of atlases: format incorrect')
                else:
                    raise ValueError('List of atlases: format incorrect')

        elif isinstance(at, dict) and ('atlas_name' in at) and \
                ('cell_types' in at):
            pass

        # Custom atlas
        elif isinstance(at, dict):
            if 'counts' not in at:
                raise ValueError('custom atlas must have a "counts" key')
            if 'number_of_cells' not in at:
                raise ValueError('custom atlas must have a "number_of_cells"'
                                 ' key')

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
            if not isinstance(at['number_of_cells'], pd.Series):
                raise ValueError('atlas["number_of_cells"] must be a series')
            if at['counts'].shape[1] != at['number_of_cells'].shape[0]:
                raise ValueError(
                    'atlas counts and number_of_cells must have the same cells')
            if (at['counts'].columns != at['number_of_cells'].index).any():
                raise ValueError(
                    'atlas counts and number_of_cells must have the same cells')
        else:
            raise ValueError(
                    'atlas must be a str, list of str, list of dict, or dict')

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
        at = self.atlas

        if isinstance(at, str):
            self.atlas = AtlasFetcher().fetch_atlas(
                    at,
                    kind='average',
                    )

        elif isinstance(self.atlas, list) or isinstance(self.atlas, tuple):
            self.atlas = AtlasFetcher().fetch_multiple_atlases(
                    at,
                    kind='average',
                    join=self.join,
                    )

        elif isinstance(at, dict) and ('atlas_name' in at) and \
                ('cell_types' in at):
            self.atlas = AtlasFetcher().fetch_atlas(
                    at['atlas_name'],
                    kind='average',
                    cell_types=at['cell_types'],
                    )

    def merge_atlas_newdata(self):
        '''Merge the averaged atlas data and the new data

        This function sets the properties:
            - n_fixed: the number of cell types in the atlas
            - n_free: the number of cell types in the new data
            - cell_types: a 1D array with the cell types of the atlas
            - features_all: a 1D array with the features that were found in
            - matrix: a 2D array with the merged counts
            - sizes: a 1D array with the sizes of each column in the matrix
            both the atlas and the new dataset.

        NOTE: is self.normalize is True, the merged count matrix is normalized
        by 1 million total counts.
        '''

        # Intersect features
        atlas_features = self.atlas['counts'].index.values
        new_data_features = self.new_data.index.values
        features = np.intersect1d(atlas_features, new_data_features)
        self.features_all = features

        # Cells types
        self.n_fixed = n_fixed = self.atlas['counts'].shape[1]
        self.n_free = n_free = self.new_data.shape[1]
        self.cell_types = self.atlas['counts'].columns.values

        # Count matrix
        L = len(features)
        N = n_fixed + n_free
        matrix = np.empty((L, N), dtype=np.float32)
        matrix[:, :n_fixed] = self.atlas['counts'].loc[features].values
        matrix[:, n_fixed:] = self.new_data.loc[features].values
        if self.normalize_counts:
            matrix *= 1e6 / matrix.sum(axis=0)
        self.matrix_all = matrix
        self.matrix = self.matrix_all

        # Cell numbers
        self.sizes = np.ones(N, np.float32)
        if self.n_cells_per_type is not None:
            self.sizes[:self.n_fixed] *= self.n_cells_per_type
        else:
            self.sizes[:self.n_fixed] = self.atlas['number_of_cells'].values.astype(np.float32)

    def select_features(self):
        '''Select features that define heterogeneity of the atlas and new data

        Returns:
            ndarray of feature names.
        '''
        # Shorten arg names
        features = self.features
        features_all = list(self.features_all)
        features_add = self.features_additional
        matrix_all = self.matrix_all

        if features is not None:
            ind_features = []
            for fea in features:
                ind_features.append(features_all.index(fea))
            self.features_selected = features
            self.matrix = matrix_all[ind_features]
            return

        n_fixed = self.n_fixed
        nf1 = self.n_features_per_cell_type
        nf2 = self.n_features_overdispersed

        ind_features = set()

        # Atlas markers
        if (nf1 > 0) and (n_fixed > 1):
            for icol in range(n_fixed):
                ge1 = matrix_all[:, icol]
                ge2 = (matrix_all[:, :n_fixed].sum(axis=1) - ge1) / (n_fixed - 1)
                fold_change = np.log2(ge1 + 0.1) - np.log2(ge2 + 0.1)
                markers = np.argpartition(fold_change, -nf1)[-nf1:]
                ind_features |= set(markers)

        # Unbiased on new data
        if nf2 > 0:
            nd_mean = matrix_all[:, n_fixed:].mean(axis=1)
            nd_var = matrix_all[:, n_fixed:].var(axis=1)
            fano = (nd_var + 1e-10) / (nd_mean + 1e-10)
            overdispersed = np.argpartition(fano, -nf2)[-nf2:]
            ind_features |= set(overdispersed)

        # Additional features
        if features_add is not None:
            tmp = pd.Series(np.arange(len(features_all)), index=features_all)
            ind_features |= set(tmp.loc[features_add].values)

        ind_features = list(ind_features)

        self.features_selected = self.features_all[ind_features]
        self.matrix = matrix_all[ind_features]

    def compute_pca(self):
        '''Compute k nearest neighbors from a matrix with fixed nodes

        Returns:
            list of lists with the first k or less indices of the neighbors for
            each free column. The length is N - n_fixed. For each now, there are
            less than k entries if no other column within the distance threshold
            were found, or if N < k.


        The algorithm proceeds as follows:
        0. take the log of the counts
        1. subtract the mean along the observation axis (N) and divide by the
        standard dev along the same axis
        2. calculate the weighted covariance matrix
        3. calculate normal PCA on that matrix
        '''
        matrix = self.matrix
        sizes = self.sizes
        n_fixed = self.n_fixed
        n_pcs = self.n_pcs

        # Test input arguments
        L, N = matrix.shape
        if len(sizes) != N:
            raise ValueError('Matrix and sizes dimensions do not match')
        if n_fixed >= N:
            raise ValueError('n_fixed larger or equal matrix number of columns')
        if n_pcs > min(L, N):
            raise ValueError('n_pcs greater than smaller matrix dimension, those eigenvalues are zero')

        # 0. take log
        matrix = np.log10(matrix + 0.1)

        # 1. standardize
        weights = 1.0 * sizes / sizes.sum()
        mean_w = matrix @ weights
        var_w = ((matrix.T - mean_w)**2).T @ weights
        std_w = np.sqrt(var_w)
        Xnorm = ((matrix.T - mean_w) / std_w).T

        # take care of non-varying components
        Xnorm[np.isnan(Xnorm)] = 0

        # 2. weighted covariance
        # This matrix has size L x L. Typically L ~ 500 << N, so the covariance
        # L x L is much smaller than N x N, hence it's fine
        cov_w = np.cov(Xnorm, fweights=sizes)

        # 3. PCA
        # rvects columns are the right singular vectors
        evals, evects = np.linalg.eig(cov_w)
        # sort by decreasing eigenvalue (explained variance) and truncate
        ind = evals.argsort()[::-1][:n_pcs]
        # NOTE: we do not actually need the eigenvalues anymore
        lvects = np.real(evects.T[ind])

        # calculate right singular vectors given the left singular vectors
        # NOTE: this is true even if we truncated the PCA via n_pcs << L
        # rvects columns are the right singular vectors
        rvects = (lvects @ Xnorm).T

        # 4. expand embedded vectors to account for sizes
        # NOTE: this could be done by carefully tracking multiplicities
        # in the neighborhood calculation, but it's not worth it: the
        # amount of overhead memory used here is small because only a few
        # principal components are used
        Ne = int(np.sum(sizes))
        rvectse = np.empty((Ne, n_pcs), np.float32)
        cell_type_expanded = []
        i = 0
        n_fixed_expanded = 0
        for isi, size in enumerate(sizes):
            if isi < n_fixed:
                cte = self.cell_types[isi]
                n_fixed_expanded += int(size)
            else:
                cte = ''
            cell_type_expanded.extend([cte] * int(size))
            for j in range(int(size)):
                rvectse[i] = rvects[isi]
                i += 1
        cell_type_expanded = np.array(cell_type_expanded)

        self.pca_data = {
            'pcs': rvects,
            'pcs_expanded': rvectse,
            'cell_type': cell_type_expanded,
            'n_atlas': n_fixed_expanded,
            }

    def compute_similarity_graph(self):
        '''Compute similarity graph from the extended PC space

        1. calculate the distance matrix by expanding atlas columns
        2. calculate neighborhoods
        3. construct similarity graph from neighborhood lists
        '''
        from scipy.spatial.distance import cdist
        import igraph as ig

        sizes = self.sizes
        n_fixed = self.n_fixed
        k = self.n_neighbors
        kout = self.n_neighbors_out_of_atlas
        metric = self.distance_metric
        threshold = self.threshold_neighborhood
        rvects = self.pca_data['pcs']
        rvectse = self.pca_data['pcs_expanded']
        Ne = len(rvectse)

        # 5. calculate distance matrix and neighbors
        # we do it row by row, it costs a bit in terms of runtime but
        # has huge savings in terms of memory since we don't need the square
        # distance matrix
        n_fixede = int(np.sum(sizes[:n_fixed]))
        neighbors = []

        # Treat things within and outside of the atlas differently
        # Atlas neighbors
        i = 0
        for isi in range(n_fixed):
            # Find the nearest neighbors in the new data
            drow = cdist(rvects[[isi]], rvects[n_fixed:], metric=metric)[0]
            ind = np.argpartition(-drow, -kout)[-kout:]

            # Discard the ones beyond threshold
            ind = ind[drow[ind] <= threshold]

            # Indices are not sorted within ind, so we need to sort them
            # in descending order by distance (more efficient in the next step)
            ind = ind[np.argsort(drow[ind])]

            for ii in range(int(sizes[isi])):
                # Internal edges
                neis = list(range(i, i+int(sizes[isi])))
                # Remove self
                neis.remove(i+ii)
                # External edges
                neis.extend(list(ind + n_fixede))
                neighbors.append(neis)
            i += int(sizes[isi])

        # New data neighbors
        for i in range(n_fixede, Ne):
            drow = cdist(rvectse[[i]], rvectse, metric=metric)[0]

            # set distance to self as a high number, to avoid self
            drow[i] = drow.max() + 1

            # Find largest k negative distances (k neighbors)
            ind = np.argpartition(-drow, -k)[-k:]

            # Discard the ones beyond threshold
            ind = ind[drow[ind] <= threshold]

            # Indices are not sorted within ind, so we need to sort them
            # in descending order by distance (more efficient in the next step)
            ind = ind[np.argsort(drow[ind])]

            neighbors.append(list(ind))

        self.neighbors = neighbors

        # Construct graph from the lists of neighbors
        edges_d = set()
        for i, neis in enumerate(neighbors):
            for n in neis:
                edges_d.add(frozenset((i, n)))

        edges = [tuple(e) for e in edges_d]
        self.graph = ig.Graph(n=Ne, edges=edges, directed=False)

    def cluster_graph(self):
        '''Compute communities from a matrix with fixed nodes

        Returns:
            None, but Averages.membership is set as an array with
            size N - n_fixed with the atlas cell types of all cells from the
            new dataset.
        '''
        import inspect
        import leidenalg

        # Check whether this version of Leiden has fixed nodes support
        opt = leidenalg.Optimiser()
        sig = inspect.getfullargspec(opt.optimise_partition)
        if 'fixed_nodes' not in sig.args:
            raise ImportError('This version of the leidenalg module does not support fixed nodes. Please update to a later (development) version')

        matrix = self.matrix
        sizes = self.sizes
        n_fixed = self.n_fixed
        clustering_metric = self.clustering_metric
        resolution_parameter = self.resolution_parameter
        g = self.graph

        L, N = matrix.shape
        n_fixede = int(np.sum(sizes[:n_fixed]))
        Ne = int(np.sum(sizes))

        # NOTE: initial membership is singletons except for atlas nodes, which
        # get the membership they have.
        initial_membership = []
        for isi in range(N):
            if isi < n_fixed:
                for ii in range(int(self.sizes[isi])):
                    initial_membership.append(isi)
            else:
                initial_membership.append(isi)

        if len(initial_membership) != Ne:
            raise ValueError('initial_membership list has wrong length!')

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

        fixed_nodes = [int(i < n_fixede) for i in range(Ne)]
        opt.optimise_partition(partition, fixed_nodes=fixed_nodes)
        membership = partition.membership[n_fixede:]

        # Convert the known cell types
        lstring = len(max(self.cell_types, key=len))
        self.membership = np.array(
                [str(x) for x in membership],
                dtype='U{:}'.format(lstring))
        for i, ct in enumerate(self.cell_types):
            self.membership[self.membership == str(i)] = ct

    def estimate_closest_atlas_cell_type(self):
        '''Estimate atlas cell type closest to each new cluster'''
        from scipy.spatial.distance import cdist

        # Use PC space
        rvectse = self.pca_data['pcs']
        n_atlas = self.pca_data['n_atlas']
        cell_types = self.pca_data['cell_type'][:n_atlas]
        L = rvectse.shape[1]

        # Extract atlas averages in PC space
        cell_types_atlas = np.unique(cell_types)
        rvectse_atlas = rvectse[:n_atlas]
        N = len(cell_types_atlas)
        avg_atl = np.empty((L, N), np.float32)
        for i, ct in enumerate(cell_types_atlas):
            # They are already replicates, take the first copy
            avg_atl[:, i] = rvectse_atlas[cell_types == ct][0]

        # Calculate averages for the new clusters
        cell_types_new = list(set(self.membership) - set(cell_types_atlas))
        rvectse_new = rvectse[n_atlas:]
        N = len(cell_types_new)
        avg_new = np.empty((L, N), np.float32)
        for i, ct in enumerate(cell_types_new):
            avg_new[:, i] = rvectse_new[self.membership == ct].mean(axis=0)

        # Calculate distance matrix between new and old in the high-dimensional
        # feature-selected space
        dmat = cdist(avg_new.T, avg_atl.T, metric='euclidean')

        # Pick the closest
        closest = np.argmin(dmat, axis=1)

        # Give it actual names
        closest = pd.Series(cell_types[closest], index=cell_types_new)

        return closest
