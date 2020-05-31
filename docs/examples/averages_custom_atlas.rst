Mapping data onto custom atlas
========================================
In this example, we will map cells onto a custom atlas, using the `Averages` class:

.. code-block:: python

   import anndata
   import northstar

   # Read in the new data to be annotated
   # Here we assume it's a loom file, but
   # of course it can be whatever format
   newdata = anndata.read_loom('...')

   # Read in the atlas with annotations
   atlas_full = anndata.read_loom('...')

   # Make sure the 'CellType' column is set
   # if it has another name, rename it
   atlas_full.obs['CellType'] = atlas_full.obs['cluster'].astype(str)

   # Subsample the atlas, we don't need
   # 1M cells to find out 5 cell types
   atlas_ave = northstar.average_atlas(
       atlas_full,
   )

   # Prepare the classifier
   # We exclude the fetal cells to focus
   # on adult tissue. To keep the fetal
   # cells, just take away the _nofetal
   # It is common to balance all cell types
   # with the same number of cells to keep
   # a high resolution in the PC space,
   model = northstar.Averages(
       atlas=atlas_ave,
       n_cells_per_type=20,
   )

   # Run the classification
   model.fit(newdata)

   # Get the inferred cell types
   cell_types_newdata = model.membership

   # Get UMAP coordinates of the atlas
   # and new data (joint embedding)
   embedding = model.embed('umap')

Notice that the `n_cells_per_type=20` is an easy way to balance the importance of each cell type in the final PC space, however it is not absolutely necessary. In fact, a larger number for some cell types will increase their weight in the PC space and therefore might increase the ability of assign cells to those types. The easiest way to use unbalanced averages is to do as follows:

.. code-block:: python

   # Let's assume you are most interested in B cells and T cells
   atlas_ave.obs['NumberOfCells'] = 20
   atlas_ave.obs.loc[['B cell', 'T cell'], 'NumberOfCells'] = 100

   model = northstar.Averages(
       atlas=atlas_ave,
   )

Notice that we omitted the `n_cells_per_type` argument in this case. The rest of the code stays the same.
