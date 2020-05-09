[![Build Status](https://travis-ci.org/northstaratlas/northstar.svg?branch=master)](https://travis-ci.org/northstaratlas/northstar)
[![License: MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![ReleaseVersion](https://img.shields.io/pypi/v/northstar?color=limegreen)](https://pypi.org/project/northstar/)
[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar.svg?type=shield)](https://app.fossa.io/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar?ref=badge_shield)
[![Documentation Status](https://readthedocs.org/projects/northstar/badge/?version=latest)](https://northstar.readthedocs.io/en/latest/?badge=latest)

![Logo](https://raw.githubusercontent.com/northstaratlas/northstar/master/docs/_static/logo.png)
# northstar
Single cell type annotation guided by cell atlases, with freedom to be queer.

## Brief description
`northstar` is a Python package to identify cell types within single cell transcriptomics datasets.
northstar's superpower is that it learns from cell atlases but still allows queer cells to make their own cluster if they want to.

Also, northstar was heavily developed during [Pride Month](https://en.wikipedia.org/wiki/Gay_pride).

## Atlas resources
![Atlas averages](https://northstaratlas.github.io/atlas_landmarks/static/logo.png)

Curated averages and subsamples from several atlases: https://northstaratlas.github.io/atlas_landmarks/

If you want us to add you cell atlas, open an issue on: https://github.com/northstaratlas/atlas_landmarks/issues

## Documentation
https://northstar.readthedocs.io

## Installation
```bash
pip install northstar
```

To automatically download and use our online atlas collection at [https://northstaratlas.github.io/atlas_averages/](https://northstaratlas.github.io/atlas_averages/), you will need to call:

```bash
pip install 'northstar[atlas-fetcher]'
```

### Dependencies
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

## Usage
See the paper below or the documentation for detailed instructions and examples. The simplest way to use `northstar` is to classify a new single cell dataset using one of the available atlases, e.g. `Darmanis_2015` on brain cells:

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
model = northstar.Averages(
        atlas='Darmanis_2015',
        n_neighbors=5,
        n_pcs=10,
        )

# Run the classifier
model.fit(new_dataset)

# Get the cluster memberships for the new cells
membership = model.membership
```

## Citation
If you use this software please cite the following paper:

Fabio Zanini\*, Bojk A. Berghuis\*, Robert C. Jones, Benedetta Nicolis di Robilant, Rachel Yuan Nong, Jeffrey Norton, Michael F. Clarke, Stephen R. Quake. **Northstar enables automatic classification of known and novel cell types from tumor samples.** bioRxiv 820928; doi: https://doi.org/10.1101/820928 

## License
`northstar` is released under the MIT license.

NOTE: The module leidenalg to perform graph-based clstering is released
under the GLP3 license. You agree with those licensing terms if you use
leidenalg within northstar.


[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2Fiosonofabio%2Fnorthstar?ref=badge_large)
