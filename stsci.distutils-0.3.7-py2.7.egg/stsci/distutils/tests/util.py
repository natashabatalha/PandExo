from __future__ import with_statement

import contextlib
import os
import shutil
import stat

try:
    reload = reload
except NameError:
    from imp import reload

from ConfigParser import ConfigParser
from distutils.ccompiler import new_compiler
from distutils.msvccompiler import MSVCCompiler
from distutils.sysconfig import customize_compiler


@contextlib.contextmanager
def open_config(filename):
    cfg = ConfigParser()
    cfg.read(filename)
    yield cfg
    with open(filename, 'w') as fp:
        cfg.write(fp)


def get_compiler_command():
    """
    Returns the name of the executable used by the default compiler on the
    system used by distutils to build C extensions.
    """

    compiler = new_compiler()
    customize_compiler(compiler)
    if isinstance(compiler, MSVCCompiler):
        compiler.initialize()
        # Normally the compiler path will be quoted as it contains spaces
        return '"%s"' % compiler.cc
    return compiler.compiler[0]


def rmtree(path):
    """
    shutil.rmtree() with error handler for 'access denied' from trying to
    delete read-only files.
    """

    def onerror(func, path, exc_info):
        if not os.access(path, os.W_OK):
            os.chmod(path, stat.S_IWUSR)
            func(path)
        else:
            raise

    return shutil.rmtree(path, onerror=onerror)
