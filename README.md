[![Build Status](https://travis-ci.org/iosonofabio/northstar.svg?branch=master)](https://travis-ci.org/iosonofabio/northstar)
[![License: MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![ReleaseVersion](https://img.shields.io/pypi/v/northstar.svg)](https://pypi.org/project/northstar/)
[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar.svg?type=shield)](https://app.fossa.io/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar?ref=badge_shield)
[![Documentation Status](https://readthedocs.org/projects/northstar/badge/?version=latest)](https://northstar.readthedocs.io/en/latest/?badge=latest)
<!--
[![Coverage Status](https://coveralls.io/repos/github/iosonofabio/northstar/badge.svg?branch=master)](https://coveralls.io/github/iosonofabio/northstar?branch=master)
-->

![Logo](docs/_static/logo.png)
# northstar
Cell type annotation guided by cell atlases, with freedom to be queer.

## Brief description
`northstar` is a Python package for single cell gene expression analysis. It computes cell clusters, which often represent cell types or cell states, similar to other
algorithms such as Louvain and Leiden community detection (graph-based) and DBSCAN and HDBSCAN (distance-based). Unlike those methods, however, northstar is
semi-supervised by a previously annotated cell atlas. Given cluster/cell type annotations from the training atlas and a new, unannotated dataset, northstar
performs clustering of the new dataset allowing cells to either belong to an atlas cluster or to form new clusters.

## Rationale
By combining the information provided by the atlas with the freedom to call new clusters, northstar tries to combine the best of both unsupervised and
supervised machine learning. Since tissues can be extremely heterogeneous, the freedom to discover new queer clusters at a reasonable computational cost
is the main strength of northstar.

Also, northstar was mostly developed during [Pride Month](https://en.wikipedia.org/wiki/Gay_pride), so we couldn't abstain from showing our support.


## Installation
For now, you can use the development version.

### Dependencies
- `numpy`
- `scipy`
- `pandas`
- `scikit-learn`
- `igraph` and `python-igraph`: this is best done by installing directly `python-igraph` via pip. That will also install the C core `iGraph` library. If you are on Windows, use the binaries as suggested on the `python-igraph` GitHub page.
- `leidenalg`: you need the develop git branch (instruction here below).

To automatically download and use our curated online atlas collection at https://iosonofabio.github.io/atlas_averages/, you will also need:
- `requests`
- `loompy`

To install `leidenalg in its development branch`:
```bash
git clone --branch develop --single-branch https://github.com/vtraag/leidenalg.git
cd leidenalg
python setup.py install
```

Once all dependencies are installed, clone this repo:
```bash
git clone https://github.com/iosonofabio/northstar.git
```
Then `cd` into it and run the setup the usual Python way:
```bash
cd northstar
python setup.py install
```

## Usage
```python
import northstar

# Choose an atlas
atlas_name = 'Darmanis_2015'

# Get a gene expression matrix of the new dataset (here a
# random matrix for simplicity)
N = 200
L = 50
new_dataset = pd.DataFrame(
    data=np.random.rand(L, N).astype(np.float32),
    index=<gene_list>,
    columns=['cell_'+str(i+1) for i in range(N)],
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
```

## Roadmap
We are planning to release on Pypi and write up a paper
to describe it.

## License notes
NOTE: The module leidenalg to perform graph-based clstering is released
under the GLP3 license. You agree with those licensing terms if you use
leidenalg within northstar.


[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar?ref=badge_large)
