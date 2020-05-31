Using Averages with a precomputed atlas
========================================
In this example, we will map cells from human brain tumors (glioblastoma) onto an existing brain atlas in our precompiled `atlas landmarks <https://northstaratlas.github.io/atlas_landmarks/>`_, using the `Averages` class:

.. code-block:: python

   import anndata
   import northstar

   # Read in the GBM data to be annotated
   # Here we assume it's a loom file, but
   # of course it can be whatever format
   newdata = anndata.read_loom('...')

   # Prepare the classifier
   # We exclude the fetal cells to focus
   # on adult tissue. To keep the fetal
   # cells, just take away the _nofetal
   model = northstar.Averages(
       atlas='Darmanis_2015_nofetal',
   )

   # Run the classification
   model.fit(newdata)

   # Get the inferred cell types
   cell_types_newdata = model.membership

   # Get UMAP coordinates of the atlas
   # and new data (joint embedding)
   embedding = model.embed('umap')
