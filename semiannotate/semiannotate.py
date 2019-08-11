# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    SemiAnnotate with atlas averages
__all__ = ['SemiAnnotate']


import numpy as np


class SemiAnnotate(object):
    '''Annotate new cell types using averages of an atlas'''

    def __init__(
        self,
        matrix,
        sizes,
        n_fixed,
        n_features_per_cell_type=30,
        n_features_overdispersed=500,
        n_pcs=20,
        n_neighbors=10,
        distance_metric='correlation',
        threshold_neighborhood=0.8,
        clustering_metric='cpm',
        resolution_parameter=0.001,
        ):
        '''Prepare the model for cell annotation

        Args:
            matrix (L x N float ndarray): matrix of data. Rows are variables
            (features/genes) and columns are observations (samples/cells).

            sizes (N array of ints): number of observations (samples/cells)
            compressed in each column.

            n_fixed (int): number of columns that are fixed. These are the first
            columns of the matrix. No knn for those columns is calculated.

            n_features_per_cell_type (int): number of features marking each fixed
            column (atlas cell type).

            n_features_overdispersed (int): number of unbiased, overdispersed
            features from the last N - n_fixed columns.

            n_pcs (int): number of principal components to keep in the weighted PCA

            n_neighbors (int): number of neighbors in the similarity graph

            distance_metric (str or function): metric to use as distance. If a
            string, it should be a metric accepted by scipy.spatial.distance.cdist.
            If a function, it should accept a (M x N)  and a (M x N1) data matrices
            as input and return a (N x N1) distance matrix. N includes both the
            fixed and the free columns, whereas N1 = N - n_fixed only includes the
            free columns.

            threshold_neighborhood (float): do not consider distances larger than this as
            neighbors

            clustering_metric (str): 'cpm' (default, Cell Potts Model) or
            'modularity'. Sets the type of partition used in the clustering
            step.

            resolution_parameter (float): number between 0 and 1 that sets
            how easy it is for the clustering algorithm to make new clusters
        '''

        self.matrix = matrix
        self.sizes = sizes
        self.n_fixed = n_fixed
        self.n_features_per_cell_type = n_features_per_cell_type
        self.n_features_overdispersed = n_features_overdispersed
        self.n_pcs = n_pcs
        self.n_neighbors = n_neighbors
        self.distance_metric = distance_metric
        self.threshold_neighborhood = threshold_neighborhood
        self.clustering_metric = clustering_metric
        self.resolution_parameter = resolution_parameter

    def _check_init_arguments(self):
        if self.n_neighbors >= self.matrix.shape[1] - 1:
            raise ValueError('n_neighbors is too high so the similarity graph is fully connected')

    def select_features(self):
        '''Select features that define heterogeneity of the atlas and new data

        Returns:
            ndarray of feature names.
        '''
        # Shorten arg names
        matrix = self.matrix
        n_fixed = self.n_fixed
        nf1 = self.n_features_per_cell_type
        nf2 = self.n_features_overdispersed

        features = set()

        # Atlas markers
        if n_fixed > 1:
            for icol in range(n_fixed):
                ge1 = matrix[:, icol]
                ge2 = (matrix[:, :n_fixed].sum(axis=1) - ge1) / (n_fixed - 1)
                fold_change = np.log2(ge1 + 0.1) - np.log2(ge2 + 0.1)
                markers = np.argpartition(fold_change, -nf1)[-nf1:]
                features |= set(markers)

        # Unbiased on new data
        nd_mean = matrix[:, n_fixed:].mean(axis=1)
        nd_var = matrix[:, n_fixed:].var(axis=1)
        fano = (nd_var + 1e-10) / (nd_mean + 1e-10)
        overdispersed = np.argpartition(fano, -nf2)[-nf2:]
        features |= set(overdispersed)
        features = list(features)

        self.features_selected = features

        self.matrix = self.matrix[features]

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

        # expand matrix to make space for duplicate_columns
        N = int(sum(self.sizes))
        matrix = np.zeros((L, N), dtype=self.matrix.dtype)
        atlas_annotations = []
        i = 0
        for isi in range(n_fixed):
            for ii in range(int(self.sizes[isi])):
                matrix[:, i] = self.matrix[:, isi]
                atlas_annotations.append(isi)
                i += 1
        matrix[:, i:] = self.matrix[:, n_fixed:]
        self.matrix_expanded = matrix
        self.atlas_annotations = atlas_annotations
        self.atlas_annotations_unique = np.unique(atlas_annotations)

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
            ind = ind[np.argsort(drow[ind])][::-1]

            nbi.extend(list(ind))

        self.neighbors = neighbors

    def compute_communities(self):
        '''Compute communities from a matrix with fixed nodes

        Returns:
            None, but SemiAnnotate.membership is set as an array of int with
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

        matrix = self.matrix_expanded
        aa = self.atlas_annotations
        aau = self.atlas_annotations_unique
        n_fixed = len(aa)
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
        tmp = set(aau)
        initial_membership = list(aa)
        i = 0
        for j in range(N - n_fixed):
            while i in tmp:
                i += 1
            initial_membership.append(i)
            tmp.add(i)
            i += 1
        del tmp

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
        membership = partition.membership

        self.membership = membership[n_fixed:]

    def compute_neighbors_unsafe(self):
        '''Compute k nearest neighbors from a matrix with fixed nodes

        Returns:
            list of lists with the first k or less indices of the neighbors for
            each free column. The length is N - n_fixed. For each now, there are
            less than k entries if no other column within the distance threshold
            were found, or if N < k.


        The algorithm proceeds as follows:
        1. whiten the matrix using weights (sizes), i.e. subtract the mean along
        the observation axis (N) and divide by the standard dev along the same axis
        2. calculate the weighted covariance matrix
        3. calculate normal PCA on that matrix
        4. calculate the distance matrix of the N1 free columns to all N columns
        5. for each free columnsm, calculate the k neighbors from the distance
        matrix, checking for threshold
        '''
        from scipy.spatial.distance import cdist

        matrix = self.matrix
        sizes = self.sizes
        n_fixed = self.n_fixed
        k = self.n_neighbors
        n_pcs = self.n_pcs
        metric = self.distance_metric
        threshold = self.threshold_neighborhood

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

        # 1. whiten
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
        lvects = evects.T[ind]

        # calculate right singular vectors given the left singular vectors
        # NOTE: this is true even if we truncated the PCA via n_pcs << L
        # rvects columns are the right singular vectors
        rvects = (lvects @ Xnorm).T

        # 4. calculate distance matrix
        # rvects is N x n_pcs. Let us calculate the end that includes only the free
        # observations and then call cdist which kills the columns. The resulting
        # matrix has dimensions N1 x N
        rvects_free = rvects[n_fixed:]
        dmat = cdist(rvects_free, rvects, metric=metric)
        dmax = dmat.max()

        # FIXME
        self.rvects_free = rvects_free

        # 5. calculate neighbors
        neighbors = [[] for n1 in range(N - n_fixed)]
        for i, drow in enumerate(dmat):
            nbi = neighbors[i]
            # set distance to self as a high number, to avoid self
            drow[i + n_fixed] = dmax + 1

            # Find largest k negative distances (k neighbors)
            ind = np.argpartition(-drow, -k)[-k:]

            # Discard the ones beyond threshold
            ind = ind[drow[ind] <= threshold]

            # Indices are not sorted within ind, so we need to sort them
            # in descending order by distance (more efficient in the next step)
            ind = ind[np.argsort(drow[ind])][::-1]

            # Take first k, keeping multiplicities in mind
            nei_missing = k
            while (nei_missing != 0) and (len(ind)):
                j = ind[-1]
                if j >= n_fixed:
                    nbi.append(j)
                    ind = ind[:-1]
                    nei_missing -= 1
                else:
                    size = int(sizes[j])
                    if size >= nei_missing:
                        nbi.extend([j] * nei_missing)
                        nei_missing = 0
                    else:
                        nbi.extend([j] * size)
                        nei_missing -= size
                        ind = ind[:-1]

        self.neighbors = neighbors

        # FIXME: just duplicate those columns to calculate the PCA
        N2 = N + 19 * n_fixed
        m2 = np.zeros((L, N2), dtype=matrix.dtype)
        for i in range(n_fixed):
            for ii in range(20):
                m2[:, 20 * i + ii] = matrix[:, i]
        m2[:, 20 * n_fixed:] = matrix[:, n_fixed:]
        from sklearn.decomposition import PCA
        Xnorm = ((m2.T - m2.mean(axis=1)) / m2.std(axis=1, ddof=0)).T
        Xnorm[np.isnan(Xnorm)] = 0
        pca = PCA(n_components=n_pcs)
        rvects2 = pca.fit_transform(Xnorm.T)
        dmat2 = cdist(rvects2, rvects2, metric=metric)
        dmax2 = dmat2.max()
        self.neighbors = []
        for i in range(20 * n_fixed, N2):
            dmat2[i, i] = dmax2 + 1
            drow = dmat2[i]
            ind = np.argpartition(-drow, -k)[-k:]
            ind = ind[drow[ind] <= threshold]
            ind = ind[np.argsort(drow[ind])][::-1]
            for ii, idx in enumerate(ind):
                if idx < 20 * n_fixed:
                    ind[ii] = idx // 20
                else:
                    ind[ii] = idx - (19 * n_fixed)
            self.neighbors.append(list(ind))
        self.rvects_free = rvects2[20 * n_fixed:]

        #import ipdb; ipdb.set_trace()

        # FIXME: add neighbors from atlas out
        self.katlas = []
        kk = 5
        for i in range(n_fixed):
            dcol = dmat[:, i]
            # Find largest k negative distances (kk neighbors) out of the atlas
            ind = np.argpartition(-dcol, -kk)[-kk:] + n_fixed
            self.katlas.append(list(ind))

    def compute_communities_unsafe(
        self,
        _self_loops=True,
        _node_sizes=True,
        ):
        '''Compute communities from a matrix with fixed nodes

        Returns:
            None, but SemiAnnotate.membership is set as an array of int with
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
        n_fixed = self.n_fixed
        clustering_metric = self.clustering_metric
        resolution_parameter = self.resolution_parameter
        neighbors = self.neighbors
        k = self.n_neighbors

        L, N = matrix.shape

        # Construct graph from weighted edges
        # NOTE: in theory we should add self-weight to the fixed columns
        from collections import Counter
        edges_d = Counter()
        for i, neis in enumerate(neighbors):
            for n in neis:
                edges_d[frozenset((i + n_fixed, n))] += 1

        # FIXME: add edges from the atlas out
        for i, neis in enumerate(self.katlas):
            for n in neis:
                edges_d[frozenset((i, n))] += self.sizes[i]

        # Self loops for atlas nodes, considered maximally connected knns
        # This is slightly less than kn because mutual nearest neighbors
        # recycle edges, so the actual expectation is kn (1 - k/2n).
        if _self_loops:
            for i in range(n_fixed):
                edges_d[(i, i)] = k * self.sizes[i] * (1 - 0.5 * k / self.sizes[i])

        edges = []
        weights = []
        for e, w in edges_d.items():
            edges.append(tuple(e))
            weights.append(w)
        g = ig.Graph(n=N, edges=edges, directed=False)

        # Edge weights
        g.es['weight'] = weights

        # FIXME
        self.edges = edges
        self.edge_weights = weights

        # Node weights
        if _node_sizes:
            g.vs['node_size'] = [int(s) for s in self.sizes]
        else:
            g.vs['node_size'] = [1] * N

        # Compute communities with semi-supervised Leiden
        # NOTE: initial membership is singletons. For fixed colunms, that is fine
        # because they are already fixed. For free columns, let them choose during
        # the clustering itself.
        fixed_nodes = [True if i < n_fixed else False for i in range(N)]
        if clustering_metric == 'cpm':
            partition = leidenalg.CPMVertexPartition(
                    g,
                    resolution_parameter=resolution_parameter,
                    node_sizes='node_size',
                    )
        elif clustering_metric == 'modularity':
            partition = leidenalg.ModularityVertexPartition(
                    g,
                    resolution_parameter=resolution_parameter,
                    node_sizes='node_size',
                    )
        else:
            raise ValueError(
                'clustering_metric not understood: {:}'.format(clustering_metric))

        opt.optimise_partition(partition, fixed_nodes=fixed_nodes)
        membership = partition.membership

        self.membership = membership[n_fixed:]

    def __call__(
            self,
            select_features=True,
            ):
        '''Run SemiAnnotate with averages of the atlas

        Args:
            select_features (bool): Whether to select features or to use the
            full data matrix. The latter is useful if a different feature
            selection was performed outside of SemiAnnotate.

        Returns:
            None, but this instance of SemiAnnotate acquired the property
            `membership` containing the cluster memberships (cell types) of the
            columns except the first n_fixed. The first n_fixed columns are
            assumes to have distinct memberships in the range [0, n_fixed - 1].
        '''
        self._check_init_arguments()

        # STEP 1
        if select_features:
            self.select_features()
        self.compute_neighbors()

        # STEP 2
        self.compute_communities()
