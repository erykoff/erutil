"""
Package:
    erutil

Sub-packages and modules:
    match:
       Ra/Dec matching tools (easy but slow)

"""
from .version import __version__

from . import match
from . import header
import des
import healpix
#from . import gaussfitter
