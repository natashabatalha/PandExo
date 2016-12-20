"""Functions for getting and saving SVN info for distribution."""


from __future__ import with_statement

import os
import subprocess

from .astutils import ImportVisitor, walk


def get_svn_version(path='.'):
    """Uses ``svnversion`` to get just the latest revision at the given
    path.
    """

    try:
        pipe = subprocess.Popen(['svnversion', path], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except OSError:
        return None

    if pipe.wait() != 0:
        return None

    return pipe.stdout.read().decode('latin1').strip()


def get_svn_info(path='.'):
    """Uses ``svn info`` to get the full information about the working copy at
    the given path.
    """

    path = os.path.abspath(path)

    try:
        pipe = subprocess.Popen(['svn', 'info', path], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        # stderr is redirected in order to squelch it.  Later Python versions
        # have subprocess.DEVNULL for this purpose, but it's not available in
        # 2.5
    except OSError:
        return 'unknown'

    if pipe.wait() != 0:
        return 'unknown'

    lines = []
    for line in pipe.stdout.readlines():
        line = line.decode('latin1').strip()
        if not line:
            continue
        if line.startswith('Path:'):
            line = 'Path: %s' % os.path.basename(path)
        lines.append(line)

    if not lines:
        return 'unknown'

    return '\n'.join(lines)
