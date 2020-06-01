.. northstar documentation master file
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: _static/logo.png
   :width: 549
   :alt: northstar logo

northstar
=================
Single cell type annotation guided by cell atlases, with freedom to be queer.


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

  pip install 'northstar[atlas-fetcher]'

Dependencies
````````````````````````
- `numpy`
- `scipy`
- `pandas`
- `scikit-learn`
- `anndata`
- `python-igraph>=0.8.0`
- `leidenalg>=0.8.0`

Optional deps to use our online atlases:

- `requests`
- `loompy`

It is recommended that you install python-igraph and leidenalg using `pip`. However, any installation (e.g. conda) that includes recent enough versions of both packages will work.

Usage example
-------------
Also see our :doc:`tutorial`.

.. code-block:: python

  import anndata
  import northstar
  
  # Choose an atlas
  atlas_name = 'Darmanis_2015'
  
  # Get an AnnData object with the new data to be annotated
  new_dataset = anndata.read_loom('...')
  # or any other format
  
  # Initialize northstar classes
  model = northstar.Averages(
          atlas=atlas,
          )
  
  # Run the classifier
  model.fit(new_dataset)
  
  # Get the cluster memberships for the new cells
  membership = model.membership


Citation
-----------------------------------
If you use this software please cite the following paper:

Fabio Zanini*, Bojk A. Berghuis*, Robert C. Jones, Benedetta Nicolis di Robilant, Rachel Yuan Nong, Jeffrey Norton, Michael F. Clarke, Stephen R. Quake. **Northstar enables automatic classification of known and novel cell types from tumor samples.** bioRxiv 820928; doi: https://doi.org/10.1101/820928 


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
   examples
   api


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
