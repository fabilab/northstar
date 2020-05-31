Create an atlas landmark
========================================
In this example, we will create an atlas landmark, i.e. a compressed approximation of a cell atlas. `northstar` can either average all cells within each cell type or take a subsample that is by default balanced across cell types, independently on their original abundance:

.. code-block:: python

   import anndata
   import northstar

   # Read in the atlas with annotations
   atlas_full = anndata.read_loom('...')

   # Make sure the 'CellType' column is set
   # if it has another name, rename it
   atlas_full.obs['CellType'] = atlas_full.obs['cluster'].astype(str)

   # Create a landmark of averages
   atlas_ave = northstar.average_atlas(
       atlas_full,
   )

   # Create a landmark of subsamples
   atlas_sub = northstar.subsample_atlas(
       atlas_full,
       n_cells=20,
   )


.. note::
   Creating average landmarks does not normalize the cells to one another before averaging. This choice reproduces the null expectation of taking all sequencing reads from those cells and lumping them together. This is not wrong but might be a little naive for some situations. You are free to normalize your atlas before using `average_atlas` to weight each cell to your liking (e.g. counts per millions). `subsample_atlas` does not suffer from this ambiguity.
