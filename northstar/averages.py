# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Atlas averages
__all__ = ['Averages']

import warnings
import numpy as np
import pandas as pd
from anndata import AnnData
import leidenalg
from .fetch_atlas import AtlasFetcher
from .cluster_with_annotations import ClusterWithAnnotations



class Averages(object):
    '''Annotate new cell types using averages of an atlas'''

    def __init__(
            self,
            atlas,
            n_cells_per_type=None,
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
            atlas (str, list of str, list of dict, dict, or AnnData): cell
             atlas to use. Generally there are two kind of choices:

             The first possibility selects the corresponding cell atlas or
             atlases from northstar's online list. The names of currently
             available dataset is here:

             https://github.com/iosonofabio/atlas_averages/blob/master/table.tsv

             (check the first column for atlas names). If a list of
             str, multiple atlases will be fetched and combined. Only features
             that are in all atlases will be kept. If you use this feature, be
             careful to not mix atlases from different species. If a list of
             dict, it merges atlases as above but you can specify what cell
             types to fetch from each atlas. Each element of the list must be a
             dict with two key-value pairs: 'atlas_name' is the atlas name, and
             'cell_types' must be a list of cell types to retain. Example:
             atlas=[{'atlas_name': 'Enge_2017', 'cell_types': ['alpha']}] would
             load the atlas Enge_2017 and only retain alpha cells. You can also
             use a dict to specify a single atlas and to retain only certain cell
             types. The format is as above, e.g. to select only alpha cells
             from Enge_2017 you can use:

             atlas={'atlas_name': 'Enge_2017', 'cell_types': ['alpha']}.

            The second possibility is to use a custom atlas (e.g. some
             unpublished data). 'atlas' must be an AnnData object with cell
             type averages ("cells") as rows and genes as columns and one cell
             metadata column 'NumberOfCells' describing the number of cells
             for each cell type. In other words:

                 adata.obs['NumberOfCells']

             must exist, and:

                adata.obs_names

            must contain the known cell types.

            n_cells_per_type (None or int): if None, use the number of cells
             per type from the atlas. Else, fix it to this number for all types.

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
        self.compute_feature_intersection()
        self._check_feature_intersection()
        self.prepare_feature_selection()
        self.select_features()
        self._check_feature_selection()
        self.merge_atlas_newdata()
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
        # Check atlas
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

        elif isinstance(at, AnnData):
            if 'NumberOfCells' not in at.obs:
                raise AttributeError(
                        'atlas must have a "NumberOfCells" obs column')

        else:
            raise ValueError('Atlas not formatted correctly')


        # Convert new data to anndata if needed
        nd = self.new_data
        if isinstance(nd, AnnData):
            pass
        elif isinstance(nd, pd.DataFrame):
            # AnnData uses features as columns, so transpose and convert
            # (the assumption is that in the DataFrame convention, rows are
            # features)
            nd = AnnData(
                X=nd.values.T,
                obs={'CellID': nd.columns.values},
                var={'GeneName': nd.index.values},
                )
            nd.obs_names = nd.obs['CellID']
            nd.var_names = nd.var['GeneName']
            self.new_data = nd
        else:
            raise ValueError(
                'New data must be an AnnData object or pd.DataFrame',
                )

        # New data could be too small to do PCA
        n_newcells, n_newgenes = self.new_data.shape
        if n_newgenes < self.n_pcs:
            warnings.warn(
                ('The number of features in the new data is lenn than ' +
                 'the number of PCs, so northstar might give inaccurate ' +
                 'results'))

        if n_newcells < self.n_pcs:
            warnings.warn(
                ('The number of cells in the new data is lenn than ' +
                 'the number of PCs, so northstar might give inaccurate ' +
                 'results'))

        if min(n_newgenes, n_newcells) < self.n_pcs:
            warnings.warn('Reducing the number of PCs to {:}'.format(
                min(n_newgenes, n_newcells)))
            self.n_pcs = min(n_newgenes, n_newcells)

        # New data could be too small for knn
        if n_newcells < self.n_neighbors + 1:
            warnings.warn(
                ('The number of cells in the new data is less than the ' +
                 'number of neighbors requested for the knn: reducing the ' +
                 'number of graph neighbors to {:}'.format(
                     max(1, n_newcells - 1)),
                 ))
            self.n_neighbors = max(1, n_newcells - 1)

        nf1 = self.n_features_per_cell_type
        nf2 = self.n_features_overdispersed
        nf3 = self.features_additional
        if not isinstance(nf1, int):
            raise ValueError('n_features_per_cell_type must be an int >= 0')
        if not isinstance(nf1, int):
            raise ValueError('n_features_overdispersed must be an int >= 0')
        if (nf1 < 1) and (nf2 < 1) and (nf3 < 1):
            raise ValueError('No features selected')

    def _check_feature_intersection(self):
        L = len(self.features_ovl)
        if L == 0:
            raise ValueError(
                ('No overlapping features in atlas and new data, are gene ' +
                 'names correct for this species?'))
        if L < 50:
            warnings.warn(
                ('Only {:} overlapping features found in atlas and new ' +
                 'data'.format(L)))

    def _check_feature_selection(self):
        L = len(self.features)
        if L == 0:
            raise ValueError(
                ('No features survived selection, check nortstar parameters'))
        if L < self.n_pcs:
            warnings.warn(
                ('Only {0} features selected, reducing PCA to {0} components'.format(L)))
            self.n_pcs = L

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

    def compute_feature_intersection(self):
        '''Calculate the intersection of features between atlas and new data'''
        # Intersect features
        self.features_atlas = self.atlas.var_names.values
        self.features_newdata = self.new_data.var_names.values
        self.features_ovl = np.intersect1d(
                self.features_atlas,
                self.features_newdata,
                )

    def prepare_feature_selection(self):
        # Cell names and types
        self.cell_types_atlas = self.atlas.obs_names
        self.cell_names_atlas = self.atlas.obs_names
        self.cell_names_newdata = self.new_data.obs_names
        ctypes_ext = []
        cnames_ext = []
        if self.n_cells_per_type is None:
            ncells_per_ct = self.atlas.obs['NumberOfCells'].astype(np.int64)
        else:
            ncells_per_ct = [self.n_cells_per_type] * self.atlas.shape[0]
        for i, ni in enumerate(ncells_per_ct):
            for ii in range(ni):
                ctypes_ext.append(self.cell_types_atlas[i])
                cnames_ext.append(self.cell_types_atlas[i]+'_{:}'.format(ii+1))
        self.cell_types_atlas_extended = ctypes_ext
        self.cell_names_atlas_extended = cnames_ext

        # Numbers
        self.n_atlas = self.atlas.shape[0]
        self.n_newdata = self.new_data.shape[0]
        self.n_total = self.n_atlas + self.n_newdata
        self.n_atlas_extended = len(self.cell_names_atlas_extended)
        self.n_total_extended = self.n_atlas_extended + self.n_newdata

        # Cell numbers
        self.sizes = np.ones(self.n_total, np.float32)
        if self.n_cells_per_type is not None:
            self.sizes[:self.n_atlas] *= self.n_cells_per_type
        else:
            self.sizes[:self.n_atlas] = self.atlas.obs['NumberOfCells'].astype(
                    np.float32)

    def select_features(self):
        '''Select features among the overlap of atlas and new data

        Returns:
            ndarray of feature names.

        '''
        features_atlas = self.features_atlas
        features_newdata = self.features_newdata
        features_ovl = list(self.features_ovl)
        features_add = self.features_additional

        n_atlas = self.n_atlas
        nf1 = self.n_features_per_cell_type
        nf2 = self.n_features_overdispersed

        features = set()

        # Atlas markers
        if (nf1 > 0) and (n_atlas > 1):
            matrix = self.atlas.X
            for icol in range(n_atlas):
                ge1 = matrix[icol]
                ge2 = (matrix.sum(axis=0) - ge1) / (n_atlas - 1)
                fold_change = np.log2(ge1 + 0.1) - np.log2(ge2 + 0.1)
                tmp = np.argsort(fold_change)[::-1]
                ind_markers_atlas = []
                for i in tmp:
                    if features_atlas[i] in features_ovl:
                        ind_markers_atlas.append(i)
                    if len(ind_markers_atlas) == nf1:
                        break
                # Add atlas markers
                features |= set(features_atlas[ind_markers_atlas])

        # Overdispersed features from new data
        if nf2 > 0:
            if nf2 >= len(features_ovl):
                features |= set(features_ovl)
            else:
                matrix = self.new_data.X
                nd_mean = matrix.mean(axis=0)
                nd_var = matrix.var(axis=0)
                fano = (nd_var + 1e-10) / (nd_mean + 1e-10)
                tmp = np.argsort(fano)[::-1]
                ind_ovd_newdata = []
                for i in tmp:
                    if features_newdata[i] in features_ovl:
                        ind_ovd_newdata.append(i)
                    if len(ind_ovd_newdata) == nf2:
                        break
                # Add overdispersed features
                features |= set(features_newdata[ind_ovd_newdata])

        # Additional features
        if features_add is not None:
            features |= (set(features_add) & set(features_ovl))

        self.features = np.array(list(features))

    def merge_atlas_newdata(self):
        '''Merge atlas data and the new data after feature selection

        NOTE: is self.normalize is True, the merged count matrix is normalized
        by 1 million total counts.
        '''
        features = self.features
        L = len(features)
        N1 = self.n_atlas
        N = self.n_total

        # This is the largest memory footprint of northstar
        matrix = np.empty((N, L), dtype=np.float32)

        # Find the feature indices for atlas
        ind_features_atlas = pd.Series(
            np.arange(len(self.features_atlas)),
            index=self.features_atlas,
            ).loc[features].values
        matrix[:N1] = self.atlas.X[:, ind_features_atlas]

        # Find the feature indices for new data
        ind_features_newdata = pd.Series(
            np.arange(len(self.features_newdata)),
            index=self.features_newdata,
            ).loc[features].values
        matrix[N1:] = self.new_data.X[:, ind_features_newdata]

        # The normalization function also sets pseudocounts
        if self.normalize_counts:
            matrix = 1e6 * (matrix.T / (matrix.sum(axis=1) + 0.1)).T

        self.matrix = matrix

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
        n_atlas = self.n_atlas
        n_pcs = self.n_pcs

        # Test input arguments
        N, L = matrix.shape
        if len(sizes) != N:
            raise ValueError('Matrix and sizes dimensions do not match')
        if n_atlas >= N:
            raise ValueError('n_fixed larger or equal matrix number of columns')
        if n_pcs > min(L, N):
            raise ValueError('n_pcs greater than smaller matrix dimension, those eigenvalues are zero')

        # 0. take log
        matrix = np.log10(matrix + 0.1)

        # 1. standardize
        weights = 1.0 * sizes / sizes.sum()
        mean_w = weights @ matrix
        var_w = weights @ ((matrix - mean_w)**2)
        std_w = np.sqrt(var_w)
        Xnorm = (matrix - mean_w) / std_w

        # take care of non-varying components
        Xnorm[np.isnan(Xnorm)] = 0

        # 2. weighted covariance
        # This matrix has size L x L. Typically L ~ 500 << N, so the covariance
        # L x L is much smaller than N x N
        cov_w = np.cov(Xnorm.T, fweights=sizes)

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
        rvects = (lvects @ Xnorm.T).T

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
            if isi < n_atlas:
                cte = self.cell_types_atlas[isi]
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
        n_atlas = self.n_atlas
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
        n_fixede = int(np.sum(sizes[:n_atlas]))
        neighbors = []

        # Treat things within and outside of the atlas differently
        # Atlas neighbors
        i = 0
        for isi in range(n_atlas):
            # Find the nearest neighbors in the new data
            drow = cdist(rvects[[isi]], rvects[n_atlas:], metric=metric)[0]
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
        clu = ClusterWithAnnotations(
                self.graph,
                self.cell_types_atlas_extended,
                resolution_parameter=self.resolution_parameter,
                metric=self.clustering_metric,
            )
        self.membership = clu.fit_transform()

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

    def embed(self, method='tsne', **kwargs):
        X = self.pca_data['pcs_expanded']
        index = list(self.cell_names_atlas_extended) + list(self.cell_names_newdata)

        if method == 'pca':
            emb = X[:, :2]
        elif method == 'tsne':
            from sklearn.manifold import TSNE
            kwargs['perplexity'] = kwargs.get('perplexity', 30)
            model = TSNE(
                n_components=2,
                **kwargs,
                )
            emb = model.fit_transform(X)
        elif method == 'umap':
            from umap import UMAP
            model = UMAP(
                n_components=2,
                **kwargs,
                )
            emb = model.fit_transform(X)

        res = pd.DataFrame(
            emb,
            index=index,
            columns=['Dimension 1', 'Dimension 2'],
            )
        res['CellType'] = list(self.cell_types_atlas_extended) + list(self.membership)
        res['Dataset'] = (['Atlas'] * self.n_atlas_extended) + (['New'] * self.n_newdata)

        return res
