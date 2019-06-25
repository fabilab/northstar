[![Build Status](https://travis-ci.org/iosonofabio/semiannotate.svg?branch=master)](https://travis-ci.org/iosonofabio/semiannotate)
[![Coverage Status](https://coveralls.io/repos/github/iosonofabio/semiannotate/badge.svg?branch=master)](https://coveralls.io/github/iosonofabio/semiannotate?branch=master)
[![License: MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![ReleaseVersion](https://img.shields.io/pypi/v/semiannotate.svg)](https://pypi.org/project/semiannotate/)
<!--
[![Documentation Status](https://readthedocs.org/projects/semiannotate/badge/?version=master)](https://semiannotate.readthedocs.io/en/master)
-->

# semiannotate
Atlas-based cell type annotation, with freedom to be queer.

## Brief description
`semiannotate` is a Python package for single cell gene expression analysis. It computes cell clusters, which often represent cell types or cell states, similar to other
algorithms such as Louvain and Leiden community detection (graph-based) and DBSCAN and HDBSCAN (distance-based). Unlike those methods, however, SemiAnnotate is
semi-supervised by a previously annotated cell atlas. Given cluster/cell type annotations from the training atlas and a new, unannotated dataset, SemiAnnotate
performs clustering of the new dataset allowing cells to either belong to an atlas cluster or to form new clusters. By combining the information provided by the atlas
with the freedom to call new clusters, SemiAnnotate tries to combine the best of both unsupervised and supervised machine learning.

## Installation
For now, you can use the development version.

### Installing dependencies
First, install the dependencies:
- `numpy` and `scipy`: use your favourite package manager, e.g. conda, pip.
- `iGraph` and `python-igraph`: this is best done by installing directly `python-igraph` via pip. That will also install the C core `iGraph` library. If you are on Windows, use the binaries as suggested on the `python-igraph` GitHub page.
- `leidenalg`: you need the develop git branch (instruction here below).

### Installing leidenalg develop branch
```bash
git clone --branch develop --single-branch https://github.com/vtraag/leidenalg.git
cd leidenalg
python setup.py install
```

### Installing semiannotate
Once all dependencies are installed, clone this repo:
```bash
git clone https://github.com/iosonofabio/semiannotate.git
```
Then `cd` into it and run the setup the usual Python way:
```bash
cd semiannotate
python setup.py install
```

## Usage
```python

# Get a gene expression matrix of atlas data (a random matrix
# here for simplicity). Each row is a feature/gene, each
# column is the average of a cluster or cell type in the atlas
Na = 10
L = 50
atlas_averages = np.random.rand(L, Na).astype(np.float32)

# Get the size of each cluster in the atlas
atlas_sizes = np.random.randint(100, size=Na)

# Get a gene expression matrix of the new dataset, aligned with
# the same features in the same order (here another random
# matrix for simplicity)
Nnew = 200
new_dataset = np.random.rand(L, Nnew).astype(np.float32)

# Concatenate the two datasets
matrix = np.hstrack([atlas_averages, new_dataset])
N = matrix.shape[1]
sizes = np.concatenate([atlas_sizes, np.ones(N, int)])

# Initialize SemiAnnotate class
sa = SemiAnnotate(
        matrix,
        sizes=sizes,
        n_fixed=Na,
        n_neighbors=5,
        n_pcs=20,
        distance_metric='correlation',
        threshold_neighborhood=0.8,
        )

# Run the classifier
sa()

# Get the cluster memberships for the new cells
membership = sa.membership
```

## Roadmap
We are planning to release semiannotate on Pypi and write up a paper
to describe it.

## License notes
NOTE: The module leidenalg to perform graph-based clstering is released
under the GLP3 license. You agree with those licensing terms if you use
leidenalg within SemiAnnotate.
