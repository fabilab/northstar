# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    SemiAnnotate top level module
__all__ = ['SemiAnnotate']


import numpy as np


class SemiAnnotate(object):
    '''Annotate new cell types using an atlas'''

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

        self.matrix = self.matrix[features]

    def compute_neighbors(self):
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
        from scipy.sparse.linalg import eigsh
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

        # 1. whiten
        weights = 1.0 * sizes / sizes.sum()
        mean_w = matrix @ weights
        var_w = (matrix**2) @ weights
        std_w = np.sqrt(var_w)
        matrix_w = ((matrix.T - mean_w) / std_w).T

        # take care of non-varying components
        matrix_w[np.isnan(matrix_w)] = 0

        # 2. weighted covariance
        # This matrix has size L x L. Typically L ~ 500 << N, so the covariance
        # L x L is much smaller than N x N, hence it's fine
        cov_w = matrix_w @ np.diag(weights) @ matrix_w.T

        # 3. PCA
        # lvects columns are the left singular vectors L x L (gene loadings)
        evals, lvects = eigsh(cov_w, k=n_pcs)

        # calculate the right singular vectors N x N given the left singular
        # vectors and the singular values svals = np.sqrt(evals)
        # NOTE: this is true even if we truncated the PCA via n_pcs << L
        # rvects columns are the right singular vectors
        svals = np.sqrt(evals)
        rvects = matrix_w.T @ lvects @ np.diag(1.0 / svals)

        # 4. calculate distance matrix
        # rvects is N x n_pcs. Let us calculate the end that includes only the free
        # observations and then call cdist which kills the columns. The resulting
        # matrix has dimensions N1 x N
        rvects_free = rvects[n_fixed:]
        dmat = cdist(rvects_free, rvects, metric=metric)
        dmax = dmat.max()

        # 5. calculate neighbors
        neighbors = [[] for n1 in range(N - n_fixed)]
        for i, drow in enumerate(dmat):
            nbi = neighbors[i]
            # set distance to self as a high number, to avoid self
            drow[i + n_fixed] = dmax + 1

            # we might not need it, but it's ok fast
            if k < len(drow):
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
                        size = sizes[j]
                        if size >= nei_missing:
                            nbi.extend([j] * nei_missing)
                            nei_missing = 0
                        else:
                            nbi.extend([j] * size)
                            nei_missing -= size
                            ind = ind[:-1]

        self.neighbors = neighbors

    def compute_communities(self,
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

        L, N = matrix.shape

        # Construct graph from weighted edges
        # NOTE: in theory we should add self-weight to the fixed columns
        from collections import Counter
        edges_d = Counter()
        for i, neis in enumerate(neighbors):
            for n in neis:
                edges_d[frozenset((i + n_fixed, n))] += 1

        # Self loops for atlas nodes, considered maximally connected knns
        # This is slightly less than kn because mutual nearest neighbors
        # recycle edges, so the actual expectation is kn (1 - k/2n).
        if _self_loops:
            for i in range(n_fixed):
                edges_d[(i, i)] = neighbors * self.sizes[i] * (1 - 0.5 * neighbors / self.sizes[i])

        edges = []
        weights = []
        for e, w in edges_d.items():
            edges.append(tuple(e))
            weights.append(w)
        g = ig.Graph(n=N, edges=edges, directed=False)

        # Edge wrights
        g.es['weight'] = weights

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
        '''Run SemiAnnotate

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
        # STEP 1
        if select_features:
            self.select_features()
        self.compute_neighbors()

        # STEP 2
        self.compute_communities()
