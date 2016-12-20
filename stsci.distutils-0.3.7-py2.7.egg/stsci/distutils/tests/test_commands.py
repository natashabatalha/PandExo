from __future__ import with_statement

import os
import sys

try:
    import numpy
except ImportError:
    numpy = None

from nose import SkipTest

from . import StsciDistutilsTestCase
from .util import get_compiler_command, open_config


class TestCommands(StsciDistutilsTestCase):
    def test_build_optional_ext(self):
        if numpy is None:
            raise SkipTest("numpy is required to run this test")
        # The test extension in the test package is already configured to be
        # "optional" by default--we'll do one test build to make sure that goes
        # smoothly
        compiler_cmd = get_compiler_command()

        _, _, exit_code = self.run_setup('build')

        # Make sure the build went successfully; a zero exit code should be
        # good enough for our purposes
        assert exit_code == 0

        # Now let's try breaking the build
        with open(os.path.join('src', 'testext.c'), 'a') as f:
            f.write('1/0')

        # We leave off the exit status from the compiler--in most cases it will
        # say "exit status 1" but that can't be guaranteed for all compilers
        msg = ('building optional extension "stsci.testpackage.testext" '
               'failed: command \'%s\' failed with exit status' % compiler_cmd)
        # Prior to Python 2.7, distutils.log output everything to stdout; now
        # warnings and errors are output to stderr
        if sys.version_info[:2] < (2, 7):
            stderr, _, exit_code = self.run_setup('build', '--force')
        else:
            _, stderr, exit_code = self.run_setup('build', '--force')
        assert exit_code == 0
        assert stderr.splitlines()[-1].startswith(msg)

        # Test a custom fail message
        with open_config('setup.cfg') as cfg:
            cfg.set('extension=stsci.testpackage.testext', 'fail_message',
                    'Custom fail message.')

        if sys.version_info[:2] < (2, 7):
            stderr, _, exit_code = self.run_setup('build', '--force')
        else:
            _, stderr, exit_code = self.run_setup('build', '--force')
        assert exit_code == 0
        assert stderr.splitlines()[-1] == 'Custom fail message.'

        # Finally, make sure the extension is *not* treated as optional if not
        # marked as such in the config
        with open_config('setup.cfg') as cfg:
            cfg.remove_option('extension=stsci.testpackage.testext',
                              'optional')

        # This error message comes out on stderr for all Python versions AFAICT
        msg = "error: command '%s' failed with exit status" % compiler_cmd
        _, stderr, exit_code = self.run_setup('build', '--force')
        assert exit_code != 0
        assert stderr.splitlines()[-1].startswith(msg)
