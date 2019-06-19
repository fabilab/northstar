# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    semiknn entry point for import
from .semiknn import compute_neighbors, compute_communities
from .version import version
__all__ = ['compute_neighbors', 'compute_communities', 'version']
