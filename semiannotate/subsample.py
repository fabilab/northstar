# vim: fdm=indent
# author:     Fabio Zanini
# date:       5/08/19
# content:    SemiAnnotate with atlas subsampling
__all__ = ['Subsample']


import numpy as np
import pandas as pd


class Subsample(object):
    '''Annotate new cell types using an even subsample of an atlas'''

    def __init__(
        self,
        matrix,
        atlas_annotations,
        n_fixed,
        n_features_per_cell_type=30,
        n_features_overdispersed=500,
        n_pcs=20,
        n_neighbors=10,
        distance_metric='correlation',
        threshold_neighborhood=0.8,
        clustering_metric='cpm',
        resolution_parameter=0.001
        ):
        '''Prepare the model for cell annotation

        Args:
            matrix (L x N float ndarray): matrix of data. Rows are variables
            (features/genes) and columns are observations (samples/cells).

            atlas_annotations (n_fixed array of ints): annotations of the first
            n_fixed columns of the matrix, corresponding to a subsample of
            the cell atlas.

            n_fixed (int): number of columns that are fixed. These are the first
            columns of the matrix. No knn for those columns is calculated.

            n_features_per_cell_type (int): number of features marking each fixed
            column (atlas cell type).

            n_features_overdispersed (int): number of unbiased, overdispersed
            features from the last N - n_fixed columns.

            n_pcs (int): number of principal components to keep in the PCA

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
        self.atlas_annotations = atlas_annotations
        self.cell_types = np.unique(atlas_annotations)
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
        aa = self.atlas_annotations
        aau = self.cell_types
        n_fixed = self.n_fixed
        nf1 = self.n_features_per_cell_type
        nf2 = self.n_features_overdispersed

        features = set()

        # Atlas markers
        if aau > 1:
            for au in aau:
                icol = (aa == au).nonzero()[0]
                li = len(icol)
                ge1 = matrix[:, icol].mean(axis=1)
                ge2 = (matrix[:, :n_fixed].sum(axis=1) - ge1 * li) / (n_fixed - li)
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

        matrix = self.matrix
        aa = self.atlas_annotations
        aau = self.cell_types
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
        membership = partition.membership[n_fixed:]

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

    def __call__(
            self,
            select_features=True,
            ):
        '''Run SemiAnnotate with subsamples of the atlas

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

        if select_features:
            self.select_features()
        self.compute_neighbors()

        self.compute_communities()
