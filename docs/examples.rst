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
