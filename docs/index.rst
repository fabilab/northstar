.. northstar documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/logo.png
   :width: 549
   :alt: northstar logo

northstar
=================
Cell type annotation guided by cell atlases, with freedom to be queer.


Brief description
--------------------------------------
`northstar` is a Python package to identify cell types within single cell transcriptomics datasets.
It uses one or more cell atlases as a baseline and assigns each cell of your dataset to either a known
cell type from the atlas(es) or to a novel cluster. northstar's superpower is that it learn from
big data (atlases) but still allows queer cells to make their own cluster if they want to.

northstar was heavily developed during `Pride Month <https://en.wikipedia.org/wiki/Gay_pride>`_.


Installation
--------------------------------------
`northstar` is a pure Python package, so you can install it using pip:

.. code-block:: bash

  pip install northstar

To automatically download and use our online atlas collection at `https://northstaratlas.github.io/atlas_averages/ <https://northstaratlas.github.io/atlas_averages/>`_, you will need to call:

.. code-block:: bash

  pip install northstar[atlas-fetcher]

Dependencies
````````````````````````
- `numpy`
- `scipy`
- `pandas`
- `scikit-learn`
- `python-igraph>=0.8.0`
- `leidenalg>=0.8.0`

Optional deps to use our online atlases:

- `requests`
- `loompy`

It is recommended that you install python-igraph and leidenalg using `pip`. However, any installation (e.g. conda) that includes recent enough versions of both packages will work.

Usage example
-------------

.. code-block:: python

  import northstar
  
  # Choose an atlas
  atlas_name = 'Darmanis_2015'
  
  # Get a gene expression matrix of the new dataset (here a
  # random matrix for simplicity)
  n_cells = 200
  n_genes = 50
  new_dataset = pd.DataFrame(
      data=np.random.rand(n_genes, n_cells).astype(np.float32),
      index=<gene_list>,
      columns=['cell_'+str(i+1) for i in range(n_cells)],
      )
  
  # Initialize northstar classes
  sa = northstar.Averages(
          atlas='Darmanis_2015',
          new_dataset,
          n_neighbors=5,
          n_pcs=10,
          )
  
  # Run the classifier
  sa()
  
  # Get the cluster memberships for the new cells
  membership = sa.membership


License
-----------------------------------
`northstar` is released under the MIT license.

NOTE: The module `leidenalg`, which is a dependency of `northstar`,
is released under the GLP3 license. You agree with those licensing
terms if you use leidenalg within northstar.


Contents
-------------
.. toctree::
   :maxdepth: 1
   :glob:

   tutorial
   api


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
