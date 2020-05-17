import pkg_resources

name = 'pyaxon'

try:
    __version__ = pkg_resources.get_distribution('pyaxon').version
except pkg_resources.DistributionNotFound:
    __version__ = ''
