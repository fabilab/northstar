Mapping data onto custom atlas
========================================
In this example, we will map cells onto a custom atlas.

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
   atlas_sub = northstar.subsample_atlas(
       atlas_full,
   )

   # Prepare the classifier
   # We exclude the fetal cells to focus
   # on adult tissue. To keep the fetal
   # cells, just take away the _nofetal
   model = northstar.Subsample(
       atlas=atlas_sub,
   )

   # Run the classification
   model.fit(newdata)

   # Get the inferred cell types
   cell_types_newdata = model.membership

   # Get UMAP coordinates of the atlas
   # and new data (joint embedding)
   embedding = model.embed('umap')
