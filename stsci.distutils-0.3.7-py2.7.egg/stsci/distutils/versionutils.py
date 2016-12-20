"""Utilities for dealing with package version info.

See also stsci.distutils.svnutils which specifically deals with adding SVN
info to version.py modules.
"""


from __future__ import with_statement

import datetime
import os
import subprocess

from .astutils import ImportVisitor, walk


VERSION_PY_TEMPLATE = """
\"\"\"This is an automatically generated file created by %(hook_function)s.
Do not modify this file by hand.
\"\"\"

__all__ = ['__version__', '__vdate__', '__svn_revision__', '__svn_full_info__',
           '__setup_datetime__']

import datetime

__version__ = %(version)r
__vdate__ = %(vdate)r
__svn_revision__ = %(svn_revision)r
__svn_full_info__ = %(svn_full_info)r
__setup_datetime__ = %(setup_datetime)r

# what version of stsci.distutils created this version.py
stsci_distutils_version = %(stsci_distutils_version)r

if '.dev' in __version__:
    def update_svn_info():
        \"\"\"Update the SVN info if running out of an SVN working copy.\"\"\"

        import os
        import string
        import subprocess

        global __svn_revision__
        global __svn_full_info__

        path = os.path.abspath(os.path.dirname(__file__))

        run_svnversion = True

        try:
            pipe = subprocess.Popen(['svn', 'info', path],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            stdout, _ = pipe.communicate()
            if pipe.returncode == 0:
                lines = []
                for line in stdout.splitlines():
                    line = line.decode('latin1').strip()
                    if not line:
                        continue
                    lines.append(line)

                if not lines:
                    __svn_full_info__ = ['unknown']
                else:
                    __svn_full_info__ = lines
            else:
                run_svnversion = False
        except OSError:
            run_svnversion = False

        if run_svnversion:
            # If updating the __svn_full_info__ succeeded then use its output
            # to find the base of the working copy and use svnversion to get
            # the svn revision.
            for line in __svn_full_info__:
                if line.startswith('Working Copy Root Path'):
                    path = line.split(':', 1)[1].strip()
                    break

            try:
                pipe = subprocess.Popen(['svnversion', path],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
                stdout, _ = pipe.communicate()
                if pipe.returncode == 0:
                    stdout = stdout.decode('latin1').strip()
                    if stdout and stdout[0] in string.digits:
                        __svn_revision__ = stdout
            except OSError:
                pass

        # Convert __svn_full_info__ back to a string
        if isinstance(__svn_full_info__, list):
            __svn_full_info__ = '\\n'.join(__svn_full_info__)


    update_svn_info()
    del update_svn_info
"""[1:]


def package_uses_version_py(package_root, package, module_name='version'):
    """Determines whether or not a version.py module should exist in the given
    package.  Returns the full path to the version.py module, regardless of
    whether it already exists.  Otherwise returns False.

    This works by checking whether or not there are any imports from the
    'version' module in the package's ``__init__.py``.

    You should write this in your package's ``__init__.py``::

        from .version import *

    """

    pdir = os.path.join(package_root, *(package.split('.')))
    init = os.path.join(pdir, '__init__.py')
    if not os.path.exists(init):
        raise Exception('Not a valid package - no __init__.py')

    try:
        visitor = ImportVisitor()
        walk(init, visitor)
    except:
        raise SyntaxError('Not able to parse %s' % init)

    found = False
    # Check the import statements parsed from the file for an import of or
    # from the svninfo module in this package
    for imp in visitor.imports:
        if imp[0] in (module_name, '.'.join((package, module_name))):
            found = True
            break
    for imp in visitor.importfroms:
        mod = imp[0]
        name = imp[1]
        if (mod in (module_name, '.'.join((package, module_name))) or
            (mod == package and name == module_name)):
            found = True
            break

    if not found:
        return False

    return os.path.join(pdir, module_name + '.py')


def clean_version_py(package_dir, package):
    """Removes the generated version.py module from a package, but only if
    we're in an SVN working copy.
    """

    pdir = os.path.join(package_root, *(package.split('.')))
    version_py = os.path.join(pdir, 'version.py')
    if not os.path.exists(svninfo):
        return

    try:
        pipe = subprocess.Popen(['svn', 'status', svninfo],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except OSError:
        return

    if pipe.wait() != 0:
        return

    # TODO: Maybe don't assume ASCII here.  Find out the best way to handle
    # this.
    if not pipe.stdout.read().decode('latin1').startswith('?'):
        return

    os.remove(version_py)


def update_setup_datetime(filename='version.py'):
    """Update the version.py with the last time a setup command was run."""

    if not os.path.exists(filename):
        return

    d = datetime.datetime.now()

    lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    with open(filename, 'w') as f:
        for line in lines:
            if not line.startswith('__setup_datetime__'):
                f.write(line)
            else:
                f.write('__setup_datetime__ = %r\n' % d)
