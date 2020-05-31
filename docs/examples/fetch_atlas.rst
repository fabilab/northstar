Fetch a precompiled atlas landmark
========================================
In this example, we will map cells onto a custom atlas, using the `Averages` class:

.. code-block:: python

   import anndata
   import northstar

   # Initialize the class
   af = northstar.AtlasFetcher()

   # Get a list of the available landmarks
   landmarks = af.list_atlases()

   # Get one of them
   atlas_sub = af.fetch_atlas(
       'Darmanis_2015',
       kind='subsample',
   )

You can also fetch multiple atlases at once. They will be merged together. Because not all genes are present in all atlases, you can decide what to do for the genes that are missing from some atlases. In this example, we keep all genes and, for each atlas, we pad the missing genes with zeros.

.. code-block:: python

   import anndata
   import northstar

   # Initialize the class
   af = northstar.AtlasFetcher()

   # Get a list of the available landmarks
   landmarks = af.list_atlases()

   # Get two atlases (merged)
   atlas_sub = af.fetch_multiple_atlases(
       ['Darmanis_2015', 'Enge_2017'],
       kind='subsample',
       join='union',
   )

