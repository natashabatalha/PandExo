"""Hooks for zest.releaser specifically for STScI software"""


import glob
import os
import shutil
import sys

from ConfigParser import ConfigParser

from setuptools.dist import Distribution
from zest.releaser.utils import ask


DEFAULT_PACKAGE_INDEX_PATH = '/eng/ssb/web/download/packages'
PACKAGE_INDEX_URL = 'http://stsdas.stsci.edu/download/packages/index'


def is_stsci_project(workingdir):
    """
    Returns True if the product being released is from STScI and is using the
    d2to1 + stsci.distutils build/install platform.

    This is determined via some basic introspection of the project layout;
    namely that it contains a setup.cfg, and that the author-email value
    contains '@stsci.edu'.  It's ham-fisted but it should do for now.
    """

    setup_cfg = os.path.join(workingdir, 'setup.cfg')
    if not os.path.exists(setup_cfg):
        return False

    cfg = ConfigParser()
    cfg.read(setup_cfg)
    if cfg.has_option('metadata', 'author-email'):
        author_email = cfg.get('metadata', 'author-email')
    elif cfg.has_option('metadata', 'author_email'):
        author_email = cfg.get('metadata', 'author_email')
    else:
        author_email = ''

    return '@stsci.edu' in author_email


def fix_dev_version_template(data):
    """
    A postreleaser.before hook to change the dev_version_template from the
    annoying default of 'x.y.z.dev0' to just 'x.y.z.dev' without the 0.
    """

    if not is_stsci_project(data['workingdir']):
        return

    data['dev_version_template'] = '%(new_version)s.dev'


def fix_sdist_format(data):
    """
    Recent versions of zest.releaser have an annoyance that it creates .zip
    sdists instead of .tar.gz.  This is supposedly to work around a bug in
    Python 2.4 with .tar.gz sdists, but none of our software supports Python
    2.4 anyways.

    Unfortunately the only way to disable this behavior, for now, is with
    monkey-patching zest.releaser.
    """

    if not is_stsci_project(data['workingdir']):
        return

    from zest.releaser.release import Releaser

    def _my_sdist_options(self):
        return ''

    Releaser._sdist_options = _my_sdist_options


def add_to_stsci_package_index(data):
    """
    A releaser.after hook to copy the source distribution to STScI's local
    package index and update the index using basketweaver.
    """

    if not is_stsci_project(data['workingdir']):
        return

    if not data['tagdir'] or not os.path.exists(data['tagdir']):
        # Do nothing if a tag checkout was not performed
        return

    if not ask('Copy source package to STScI package index'):
        return

    package_path = DEFAULT_PACKAGE_INDEX_PATH
    if not os.path.exists(package_path):
        package_path = ''

    question = 'Path to package directory'
    if package_path:
        # A default exists; let the user know
        question += ' [%s]' % package_path
    question += ': '

    answer = ''
    while not answer:
        try:
            answer = raw_input(question).strip()
            if not answer:
                if package_path:
                    # The user simple pressed enter, so use the supplied
                    # default
                    answer = package_path
                else:
                    continue
            if not os.path.exists(answer):
                print ('The supplied path %s does not exist.  Please enter a '
                       'different path or press Ctrl-C to cancel.' % answer)
            if not os.access(answer, os.W_OK):
                print ('The supplied path %s is not writeable.  Either change '
                       'the permissions of the directory or have someone '
                       'grant you access and try again, enter a different '
                       'directory, or press Ctrl-C to cancel.' % answer)
            package_path = answer
            break
            # The default was not supplied, so keep asking
        except KeyboardInterrupt:
            return

    # A tag checkout was made and an sdist created, this is where it would be
    # (the sdist is a .zip on Windows, a .tar.gz elsewhere--normally this
    # should be .tar.gz; don't make releases on Windows)
    sdist_file = ''
    while not sdist_file:
        try:
            sdist_file = glob.glob(os.path.join(data['tagdir'], 'dist',
                                                '*.tar.gz'))[0]
        except IndexError:
            try:
                sdist_file = glob.glob(os.path.join(data['tagdir'], 'dist',
                                                    '*.zip'))[0]
            except IndexError:
                try:
                    print (
                        "Could not find a source distribution in %s; did you "
                        "do a source checkout for upload?  If possible, try "
                        "to cd to %s and manually create a source "
                        "distribution by running `python setup.py sdist`.  "
                        "Then press enter to try again (or hit Ctrl-C to "
                        "cancel).  Go ahead, I'll wait..." %
                        (data['tagdir'], data['tagdir']))
                    raw_input()
                except KeyboardInterrupt:
                    return

    # Almost ready go to--now we just need to check if basketweaver is
    # available, and get it if not.
    try:
        import basketweaver.makeindex
    except ImportError:
        # Use setuptools' machinery to fetch a package and add it to the path;
        # we could do this without using setuptools directly, but it would
        # basically end up reimplementing much of the same code.
        dist = Distribution({'dependency_links': [PACKAGE_INDEX_URL]})
        try:
            dist.fetch_build_eggs(['basketweaver'])
        except:
            # There are so many things that could possibly go wrong here...
            print ('Failed to get basketweaver, which is required to rebuild '
                   'the package index.  To manually complete the release, '
                   'install basketweaver manually, then copy %s into %s, cd '
                   'to %s, and then run `makeindex *`, where makeindex is the '
                   'command installed by basketweaver.' %
                   (sdist_file, package_path, package_path))
        import basketweaver.makeindex

    # Now we should have everything we need...
    shutil.copy(sdist_file, package_path)
    old_cwd = os.getcwd()
    os.chdir(package_path)
    try:
        basketweaver.makeindex.main(glob.glob('*'))
    finally:
        os.chdir(old_cwd)

    print 'Finished adding package to %s.' % PACKAGE_INDEX_URL
