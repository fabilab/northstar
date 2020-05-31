Tinkering with northstar's internals
========================================
The easiest way to understand northstar's internal code is to use `ipython` or `jupyter` notebooks, and use the amazing double question mark magic on the `fit` method:

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
   model = northstar.Subsample(
       atlas='Darmanis_2015_nofetal',
   )

   # Peek inside
   model.fit??

This will show you the order of operations. You can use the same double question mark on each of them to understand in detail what's going on.

Notice also that `northstar` is a pure Python package. This means that, provided you installed the dependencies, you can download the source code from pypi or GitHub, put into a folder on your machine, add that folder to your `PYTHONPATH`, and modify the code at will without needing to recompile anything. In fact, you are encouraged to do so and build even better classifiers!
