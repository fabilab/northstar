Tutorial
========
Thank you for looking into `northstar`! This tutorial guides you through a typical use of
the package to annotate your single cell dataset based on one or more cell atlases. At
the end of the tutorial, you should be able to navigate the other docs yourself (see api_)

Starting point: atlas landmarks
--------------------------------
To transfer cell types from an atlas, you need a few cells or averages for each cell type
within the atlas. We call these **atlas landmarks**. To keep things simple, in this tutorial
we use precomputed landmarks from our `sister project <https://iosonofabio.github.io/atlas_landmarks/>`_.
Alternatively, you can use a custom atlas with its landmarks, see the api_ documentation of
classes :class:`.Averages` and :class:`.Subsample` for that advanced usage.


Preparing your single cell dataset
--------------------------------------
TODO

.. toctree::
   :hidden:
   :maxdepth: 2
   :glob:

   examples/*


.. _api:
