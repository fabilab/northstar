Using an external tool for harmonization
========================================
In this example, we will take a neighborhood graph obtained with an external tool such as scVI and cluster the cells using northstar's atlas-aware clustering algorithm.

The key class is `ClusterWithAnnotations`, which takes two required arguments:
- the external `graph` must satisfy a specific order of cells. The first cells (starting from 0) must be the annotated atlas cells, the last ones must be the cells to be annotated
- the `annotations` array/list contains the cell types of the atlas cells. The length of this array/list will be used to infer which cells in your graph are annotated atlas cells (the first ones, from 0 to the length of this array - 1) and which ones are to be annotated (the other ones, i.e. the last ones).

.. code-block:: python

   import anndata
   import northstar

   # Read in the GBM data to be annotated
   # Here we assume it's a loom file, but
   # of course it can be whatever format
   newdata = anndata.read_loom('...')

   # The variable graph contains a sparse
   # adjacency matrix between cells. In
   # other words, graph[i, j] is nonzero
   # if cells i and j are neighbors. graph
   # was computed using an external tool
   print(graph)

   # The variable annotations contains an
   # array/list of cell types for the atlas
   # cells, which are cells 0 to n_atlas - 1
   # in the graph
   print(annotations)
   n_atlas = len(annotations)

   # Prepare the clustering class
   model = northstar.ClusterWithAnnotations(
       graph=graph,
       annotations=annotations,
       )

   # Run the classification
   model.fit(newdata)

   # Get the inferred cell types
   cell_types_newdata = model.membership

Although this example uses a sparse adjacency matrix, you can also use an `igraph` graph instead. A typical way to create a graph from a list of `edges` is:

.. code-block:: python

   import igraph as ig

   # edges is a list of pairs, e.g.
   # edges = [(0, 1), (0, 3)]
   # would indicate that cells 0 and 1
   # are neighbors, 0 and 3 are also
   # neighbors, while cell 2 has no
   # neighbors
   print(edges)

   graph = igraph.Graph(
       n=n_atlas+n_newcells,
       edges=edges,
       )

   # Prepare the clustering class
   model = northstar.ClusterWithAnnotations(
       graph=graph,
       annotations=annotations,
       ) 

The rest of the code stays the same.
