# vim: fdm=indent
# author:     Fabio Zanini
# date:       5/08/19
# content:    Cluster with annotations. This is used in Subsample and can be
#             used independently with external similarity graphs and annotations.
__all__ = ['ClusterWithAnnotations']

import numpy as np
import scipy as sp
import igraph as ig
import leidenalg


class ClusterWithAnnotations(object):
    '''Cluster a similarity graph of cells with fixed atlas annotations'''

    def __init__(
            self,
            graph,
            annotations,
            resolution_parameter=0.001,
            metric='cpm',
            ):
        '''Prepare the clustering with annotations


        Args:
            graph (igraph.Graph or adjecency matrix): the
                similarity graph connecting cells. The first N cells are
                considered part of the atlas, where N is the length of the
                annotations argument.

            annotations (list or numpy ndarray): cell type annotations for
                the atlas cells, assumed to be the first one in the graph.

            resolution_parameter (float): number between 0 and 1 that sets
             how easy it is for the clustering algorithm to make new clusters

            metric (str): 'cpm' (default, Cell Potts Model) or 'modularity'.
                Sets the type of partition used in the clustering step.

        '''
        self.graph = graph
        self.annotations = annotations
        self.resolution_parameter = resolution_parameter
        self.metric = metric

    def _parse_graph(self):
        '''Convert whatever input for the graph attribute into igraph.Graph'''
        arg = self.graph
        if isinstance(arg, ig.Graph):
            return

        if isinstance(arg, np.ndarray) and (arg.shape[0] == arg.shape[1]):
            n_nodes = arg.shape[0]
            edges = list(zip(*(arg.nonzero())))
        elif sp.sparse.issparse(arg):
            arg = arg.tocoo()
            n_nodes = arg.shape[0]
            edges = list(zip(arg.row, arg.col))
        else:
            raise ValueError('Graph format not recognized')

        self.graph = ig.Graph(n=n_nodes, edges=edges)

    def fit(self):
        '''Compute communities from a matrix with fixed nodes

        Returns:
            None, but the membership attribute is set as an array of int with
            size N - n_fixed with the community/cluster membership of all
            columns except the first n fixed ones.
        '''
        self._parse_graph()

        aa = self.annotations
        n_fixed = len(aa)
        g = self.graph
        N = g.vcount()

        opt = leidenalg.Optimiser()
        is_membership_fixed = [int(i < n_fixed) for i in range(N)]

        # NOTE: initial membership is singletons except for atlas nodes, which
        # get the membership they have.
        aau = list(np.unique(aa))
        aaun = len(aau)
        initial_membership = []
        for j in range(N):
            if j < n_fixed:
                mb = aau.index(aa[j])
            else:
                mb = aaun + (j - n_fixed)
            initial_membership.append(mb)

        if self.metric == 'cpm':
            partition = leidenalg.CPMVertexPartition(
                    g,
                    resolution_parameter=self.resolution_parameter,
                    initial_membership=initial_membership,
                    )
        elif self.metric == 'modularity':
            partition = leidenalg.ModularityVertexPartition(
                    g,
                    resolution_parameter=self.resolution_parameter,
                    initial_membership=initial_membership,
                    )
        else:
            raise ValueError(
                'clustering_metric not understood: {:}'.format(self.metric))

        # Run modified Leiden here
        opt.optimise_partition(partition,
                               is_membership_fixed=is_membership_fixed)

        # Exctract result
        membership = partition.membership[n_fixed:]

        # Convert the known cell types
        lstring = len(max(aau, key=len))
        self.membership = np.array(
                [str(x) for x in membership],
                dtype='U{:}'.format(lstring))
        for i, ct in enumerate(aau):
            self.membership[self.membership == str(i)] = ct

    def fit_transform(self):
        '''Cluster graph and return memberships'''
        self.fit()
        return self.membership
