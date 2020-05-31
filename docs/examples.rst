Examples
--------------------------------------
Northstar has two main classes to use averages or subsamples of cell atlases:

- :doc:`Averages with precomputed atlas <examples/averages_precomputed_atlas>`
- :doc:`Subsample with precomputed atlas <examples/subsample_precomputed_atlas>`

You can use a custom atlas:

- :doc:`Averages with custom atlas <examples/averages_custom_atlas>`
- :doc:`Subsample with custom atlas <examples/subsample_custom_atlas>`

You can also harmonize your atlas and target dataset (to be annotated) with another tool and then use northstar for clustering only. The advantage is that northstar's clustering algorithm is aware of the atlas annotations, therefore it is guaranteed to neither split not merge atlas cell types:

- :doc:`External data harmonization <examples/external_harmonization>`

You can also use northstar just as an API interface to our precompiled list of annotated atlases. This can be used to download averages and subsamples (we call them **atlas landmarks**) and use them to do whatever you want (e.g. classify using another tool, harmonize, look up marker genes, etc):

- :doc:`Fetch a precompiled atlas landmark <examples/fetch_atlas>`

Another use of northstar is for its ability to compress large datasets into cell type averages or subsamples. We call this operation **creating an atlas landmark**:

- :doc:`Compress an atlas <examples/compress_atlas>`

Finally, you might want to play with northstar's internals and take inspiration to build your very own classifer, data harmonizer, clustering algorithm, feature selection tool, or whatever else. This is probably only interesting for **advanced users**:

- :doc:`Tinkering with northstar's internals <examples/internals>`

