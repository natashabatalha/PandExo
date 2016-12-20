try:
    from .version import (__version__, __svn_revision__, __svn_full_info__,
                          __setup_datetime__)
except ImportError:
    # If this happens, I hope this import is happening in the
    # stsci.distutils setup.py because we are about to use
    # stsci.distutils.hooks to install stsci.distutils
    __version__ = 'unspecified during initial install'
    __svn_revision__ = ''
    __svn_full_info__ = ''
    __setup_datetime__ = None
