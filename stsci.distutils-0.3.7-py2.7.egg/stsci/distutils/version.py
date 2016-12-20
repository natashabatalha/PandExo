"""This is an automatically generated file created by stsci.distutils.hooks.version_setup_hook.
Do not modify this file by hand.
"""

__all__ = ['__version__', '__vdate__', '__svn_revision__', '__svn_full_info__',
           '__setup_datetime__']

import datetime

__version__ = '0.3.7'
__vdate__ = 'unspecified'
__svn_revision__ = '28700'
__svn_full_info__ = 'Path: stsci.distutils-0.3.7-eLbs04\nWorking Copy Root Path: /tmp/stsci.distutils-0.3.7-eLbs04\nURL: https://svn.stsci.edu/svn/ssb/stsci_python/stsci.distutils/tags/0.3.7\nRepository Root: https://svn.stsci.edu/svn/ssb/stsci_python\nRepository UUID: fe389314-cf27-0410-b35b-8c050e845b92\nRevision: 28700\nNode Kind: directory\nSchedule: normal\nLast Changed Author: embray\nLast Changed Rev: 28700\nLast Changed Date: 2013-12-23 18:22:33 -0500 (Mon, 23 Dec 2013)'
__setup_datetime__ = datetime.datetime(2016, 12, 20, 10, 56, 37, 848063)

# what version of stsci.distutils created this version.py
stsci_distutils_version = '0.3.7'

if '.dev' in __version__:
    def update_svn_info():
        """Update the SVN info if running out of an SVN working copy."""

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
            __svn_full_info__ = '\n'.join(__svn_full_info__)


    update_svn_info()
    del update_svn_info
