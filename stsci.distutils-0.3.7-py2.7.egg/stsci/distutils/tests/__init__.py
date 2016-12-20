import os
import shutil
import subprocess
import sys
import tempfile

import nose

from .util import reload, rmtree


TESTPACKAGE_URL = ('https://svn.stsci.edu/svn/ssb/stsci_python/'
                   'stsci.distutils/trunk/stsci/distutils/tests/testpackage')
TESTPACKAGE_REV = '17597'  # The last known 'good' revision of this package


class StsciDistutilsTestCase(object):
    @classmethod
    def setup_class(cls):
        cls.wc_dir = tempfile.mkdtemp(prefix='stsci-distutils-test-')
        try:
            p = subprocess.Popen(['svn', '-r', TESTPACKAGE_REV,
                                  '--non-interactive', '--trust-server-cert',
                                  'checkout', TESTPACKAGE_URL, cls.wc_dir],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
        except OSError, e:
            raise nose.SkipTest('svn unavailable to checkout out test '
                                'package: %s' % e)

        if p.wait() != 0:
            raise nose.SkipTest('svn failed to check out the test package: '
                                '%s; tests will not be able to run' %
                                p.stderr.read().decode('latin1'))

    @classmethod
    def teardown_class(cls):
        rmtree(cls.wc_dir)

    def setup(self):
        self.temp_dir = tempfile.mkdtemp(prefix='stsci-distutils-test-')
        self.package_dir = os.path.join(self.temp_dir, 'testpackage')
        shutil.copytree(self.wc_dir, self.package_dir)
        self.oldcwd = os.getcwd()
        os.chdir(self.package_dir)

        # We need to manually add the test package's path to the stsci
        # package's __path__ since it's already been imported.
        if 'stsci' in sys.modules:
            # Clean the existing __path__ up
            reload(sys.modules['stsci'])
            sys.modules['stsci'].__path__.insert(
                0, os.path.join(self.package_dir, 'stsci'))

    def teardown(self):
        os.chdir(self.oldcwd)
        # Remove stsci.testpackage from sys.modules so that it can be freshly
        # re-imported by the next test
        for k in list(sys.modules):
            if k == 'stsci.testpackage' or k.startswith('stsci.testpackage.'):
                del sys.modules[k]
        rmtree(self.temp_dir)

    def run_setup(self, *args):
        return self._run_cmd(sys.executable, ('setup.py',) + args)

    def run_svn(self, *args):
        return self._run_cmd('svn', args)

    def _run_cmd(self, cmd, args):
        """
        Runs a command, with the given argument list, in the root of the test
        working copy--returns the stdout and stderr streams and the exit code
        from the subprocess.
        """

        os.chdir(self.package_dir)
        p = subprocess.Popen([cmd] + list(args), stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)

        streams = tuple(s.decode('latin1').strip() for s in p.communicate())

        return (streams) + (p.returncode,)
