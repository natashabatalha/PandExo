from __future__ import with_statement


import glob
import os
import sys
import tarfile
import time
import zipfile

from datetime import datetime
from setuptools import Distribution

try:
    import numpy
except ImportError:
    numpy = None

from nose import SkipTest

from . import StsciDistutilsTestCase, TESTPACKAGE_REV
from .util import reload, get_compiler_command, open_config, rmtree


VERSION = '0.1.dev' + TESTPACKAGE_REV


class TestHooks(StsciDistutilsTestCase):
    def test_setup_py_version(self):
        """
        Test that the `./setupy.py --version` command returns the correct
        value without balking.
        """

        self.run_setup('egg_info')
        stdout, _, _ = self.run_setup('--version')
        assert stdout == VERSION

    def test_version_with_rev(self):
        """Test that the version string contains the correct SVN revision."""

        # Build the package
        self.run_setup('build')
        self.run_setup('sdist')

        import stsci.testpackage

        assert hasattr(stsci.testpackage, '__version__')
        assert stsci.testpackage.__version__ == VERSION

        assert hasattr(stsci.testpackage, '__svn_revision__')
        assert stsci.testpackage.__svn_revision__ == TESTPACKAGE_REV

        filenames = [os.path.join('dist',
                                  'stsci.testpackage-%s.%s' % (VERSION, ext))
                     for ext in ('tar.gz', 'zip')]

        assert os.path.exists(filenames[0]) or os.path.exists(filenames[1])

    def test_release_version(self):
        """
        Ensure that the SVN revision is not appended to release versions
        (i.e. not ending with '.dev'.
        """

        with open_config('setup.cfg') as cfg:
            cfg.set('metadata', 'version', '0.1')

        self.run_setup('egg_info')
        stdout, _, _ = self.run_setup('--version')
        assert stdout == '0.1'

    def test_inline_svn_update(self):
        """Test version.py's capability of updating the SVN info at runtime."""

        self.run_setup('build')

        import stsci.testpackage

        assert hasattr(stsci.testpackage, '__svn_revision__')
        assert stsci.testpackage.__svn_revision__ == TESTPACKAGE_REV

        with open('TEST', 'w') as f:
            # Create an empty file
            pass

        self.run_svn('add', 'TEST')
        # The working copy has been modified, so now svnversion (which is used
        # to generate __svn_revision__) should be the revision + 'M'
        reload(stsci.testpackage.version)
        reload(stsci.testpackage)

        assert stsci.testpackage.__svn_revision__ == TESTPACKAGE_REV + 'M'

    def test_setup_datetime(self):
        """
        Test that the setup datetime is present, and is updated by subsequent
        setup.py runs.
        """

        # Build the package
        self.run_setup('build')

        import stsci.testpackage

        assert hasattr(stsci.testpackage, '__setup_datetime__')
        prev = stsci.testpackage.__setup_datetime__
        now = datetime.now()
        # Rebuild
        # So that there's less chance for ambiguity
        time.sleep(1)
        self.run_setup('build')

        reload(stsci.testpackage.version)
        reload(stsci.testpackage)

        import stsci.testpackage

        assert hasattr(stsci.testpackage, '__setup_datetime__')
        assert stsci.testpackage.__setup_datetime__ > now
        assert stsci.testpackage.__setup_datetime__ > prev

    def test_numpy_extension_hook(self):
        """Test basic functionality of the Numpy extension hook."""

        if numpy is None:
            raise SkipTest("numpy is required to run this test")

        compiler_cmd = get_compiler_command()

        stdout, _, _ = self.run_setup('build')
        for line in stdout.splitlines():
            # Previously this used shlex.split(), but that's broken for unicode
            # strings prior to Python 3.x, and it doesn't matter too much since
            # we only care about the first argument
            args = line.split()
            if not args:
                continue
            if args[0] != compiler_cmd:
                continue

            # The first output from the compiler should be an attempt to
            # compile a c file to an object, so that should include all the
            # include paths.  This is of course not universally true, but it
            # should hold true for this test case
            for path in [numpy.get_include()]:
            #for path in [numpy.get_include(), numpy.get_numarray_include()]:
                assert '-I' + path in args
            break

        # And for the heck of it, let's ensure that this doesn't happen if
        # 'numpy' is not listed in include_dirs
        with open_config('setup.cfg') as cfg:
            cfg.remove_option('extension=stsci.testpackage.testext',
                              'include_dirs')

        rmtree('build')

        stdout, _, _ = self.run_setup('build')
        for line in stdout.splitlines():
            args = line.split()
            if not args:
                continue
            if args[0] != compiler_cmd:
                continue
            for path in [numpy.get_include()]:
            #for path in [numpy.get_include(), numpy.get_numarray_include()]:
                assert '-I' + path not in args

    def test_glob_data_files(self):
        """
        Test the glob_data_files hook by ensuring that all the correct data
        files are included in the source distribution, and that they are
        installed to the correct location in the package.
        """

        data_files = os.path.join('stsci', 'testpackage', 'data_files')

        # First test the source distribution
        self.run_setup('sdist')

        # There can be only one
        try:
            tf = glob.glob(os.path.join('dist', '*.tar.gz'))[0]
        except IndexError:
            # No tar.gz files found?  On Windows sdist creates a .zip file, so
            # let's look for that
            tf = glob.glob(os.path.join('dist', '*.zip'))[0]
            # If that failed then I don't know--I guess the test should fail

        if tf.endswith('.tar.gz'):
            tf = tarfile.open(tf)
            # Tarfiles created by sdist kindly place all contents in a
            # top-level directory with the same name as the file minus
            # extension, so as to kindly not bomb you when you extract it.  But
            # we don't care about that top level directory
            names = ['/'.join(p.split('/')[1:]) for p in tf.getnames()]
        else:
            with zipfile.ZipFile(tf) as zf:
                names = ['/'.join(p.split('/')[1:]) for p in zf.namelist()]

        # Sdists should place the data_files at the root, just like in the
        # normal source layout; even files that aren't normally installed
        # should be included
        for filename in ['a.txt', 'b.txt', 'c.rst']:
            # Don't use os.path.join -- zipfile/tarfile always return paths
            # with / as path sep
            assert ('data_files/' + filename) in names

        # Now we test that data_files go to the right place in various install
        # schemes
        def get_install_lib(args):
            # This helper uses the distutils/setuptools machinery to determine
            # where a command will install files based on the arguments passed
            # to setup.py
            dist = Distribution({'script_args': args})
            dist.parse_command_line()
            install_cmd = dist.get_command_obj('install')
            install_cmd.ensure_finalized()
            return install_cmd.install_lib

        def test_install_scheme(args):
            if numpy is None:
                raise SkipTest("numpy is required to run this test")
            # This general code should work to test the files in a variety of
            # install schemes depending on args
            if os.path.exists('temp'):
                rmtree('temp')
            install_lib = get_install_lib(args)
            os.makedirs(install_lib)
            old_pythonpath = os.environ.get('PYTHONPATH')
            # For a setuptools/easy_install-stype install to an alternate
            # prefix we have to have the new install dir on the PYTHONPATH or
            # easy_install will balk
            os.environ['PYTHONPATH'] = (
                install_lib + os.pathsep +
                (old_pythonpath if old_pythonpath else ''))

            try:
                self.run_setup(*(args + ['--record=installed.txt']))
            finally:
                if old_pythonpath is not None:
                    os.environ['PYTHONPATH'] = old_pythonpath

            found_files = 0
            with open('installed.txt') as f:
                # installed.txt, however, contains OS-specific paths
                for line in f:
                    for name in ['a.txt', 'b.txt', 'c.rst']:
                        if line.strip().endswith(os.sep + name):
                            found_files += 1
            assert found_files == 2

        test_install_scheme(['install', '--prefix=temp'])
        test_install_scheme(['install', '--root=temp'])
        test_install_scheme(['install', '--install-lib=temp'])
