# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    northstar entry point for import
from .averages import Averages
from .subsample import Subsample
from .fetch_atlas import AtlasFetcher
from ._version import version


# Check whether this version of Leiden has fixed nodes support
import inspect as _inspect
import leidenalg as _leidenalg
_opt = _leidenalg.Optimiser()
_sig = _inspect.getfullargspec(_opt.optimise_partition)
if 'fixed_nodes' not in _sig.args:
    raise ImportError('This version of the leidenalg module does not support fixed nodes. Please update to a later (development) version')


__all__ = ['Aveages', 'Subsample', 'AtlasFetcher', 'version']
