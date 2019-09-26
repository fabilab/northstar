Tutorial
========
Thank you for looking into `northstar`! This tutorial guides you through a typical use of
the package to annotate your single cell dataset based on one or more cell atlases. At
the end of the tutorial, you should be able to navigate the other docs yourself (see
:doc:`/api`).

Start off: atlas landmarks
--------------------------------
To transfer cell types from an atlas, you need a few cells or averages for each cell type
within the atlas. We call these **atlas landmarks**. To keep things simple, in this tutorial
we use precomputed landmarks from our `sister project <https://iosonofabio.github.io/atlas_landmarks/>`_.
Alternatively, you can use a custom atlas with its landmarks, see the :doc:`api` documentation of
classes :class:`.Averages` and :class:`.Subsample` for that advanced usage.

For this tutorial, we will use the atlas `Darmanis_2015 <https://www.pnas.org/content/112/23/7285>`_.


Prepare your single cell dataset
--------------------------------------
Then we need to prepare the new single cell dataset to annotate. `northstar` accepts two formats:

1. a `pandas.DataFrame` with features/genes as rows and cells as columns.
   
2. an `anndata.Anndata` object. Turns out `Anndata <https://anndata.readthedocs.io/en/stable/>`_
   enforces the opposite convention, so the rows/observation_names must be the cells and the
   columns/variable_names must be the features/genes.
   
.. note::
   `northstar` will take the intersection of your features names and the atlas features to
   assign cell types. Most atlases use gene names instead of EnsemblIDs or other names, so
   make sure you do the same.

For this tutorial, we will use pandas, so let's load the data e.g. from a CSV file:

.. code-block:: python

  # Assume the CSV already has genes as rows
  dataset = pd.read_csv('mydataset.csv', sep=',', index_col=0)

If the dataset is in a `loom <http://loompy.org/>`_ file, we have to set the index and columns:

.. code-block:: python

  import loompy

  # Assume the gene name is in an attribute called GeneName and the cell name in CellID
  with loompy.connect('mydataset.loom') as ds:
      dataset = pd.DataFrame(
          data=ds[:, :],
          index=ds.ra['GeneName'],
          columns=ds.ca['CellID'],
          )

Create an instance of northstar
-------------------------------
Let's assume you want to use the `.Subsample` class. You can create an instance of the class
as follows:

.. code-block:: python

  import northstar

  model = northstar.Subsample(
      atlas='Darmanis_2015',
      new_data=dataset,
      )

Call the cell classifier
------------------------
This is where the actual computations happen. A `northstar` instance can be called like a
normal function:

.. code-block:: python

  model()

Extract the result
------------------
The result of the cell type assignment can be extracted by the following command:

.. code-block:: python

  cell_types = model.membership

This is a numpy array with the same length and order as your cells.

Optional: closest atlas cell type
---------------------------------
Sometimes you get novel clusters that do not match any atlas cell type. To start identifying
those clusters, you can ask `northstar` what known atlas cell type they are most similar to.
Here's the code to do that:

.. code-block:: python

  closest_cell_types = model.estimate_closest_atlas_cell_type()

The output is a `pandas.Series` with the novel clusters as index and the closest atlas cell
types as values.

Where to go from here
---------------------
This concludes this short tutorial that showcases the main usage of `northstar`. You can
browse the :doc:`/api` page for more detailed information such as how to combine atlases,
specify only certain cell types within one or multiple atlases, and use custom atlases.

We hope `northstar` helps you understand your tissue sample and do not hesitate to open an
`issue on github <https://github.com/iosonofabio/northstar/issues>`_ if you have trouble.
If `northstar` was useful for a publication, please consider citing us at `PAPER MISSING`.

.. toctree::
   :hidden:
   :maxdepth: 2
   :glob:

