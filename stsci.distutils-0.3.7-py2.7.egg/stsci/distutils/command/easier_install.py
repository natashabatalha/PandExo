import os

try:
    from ConfigParser import ConfigParser
except ImportError:
    # This is necessary for when stsci.distutils is being bootstrapped by other
    # setup.pys in stsci_python--otherwise the ConfigParser import would be
    # converted by d2to1
    from configparser import ConfigParser
from distutils import log

import pkg_resources

from setuptools.command.easy_install import easy_install
from setuptools.package_index import PackageIndex


def distro_from_setup_cfg(filename):
    """
    Read a source checkout's distutils2 setup.cfg and create a Distribution for
    that checkout.

    filename can either be the path to the setup.cfg itself, or checkout
    directory containing the setup.cfg.
    """

    if os.path.isdir(filename):
        path = filename
        filename = os.path.join(filename, 'setup.cfg')
        if not os.path.exists(filename):
            return None
    else:
        path, basename = os.path.split(filename)
        if basename != 'setup.cfg':
            return None
    cfg = ConfigParser()
    cfg.read(filename)
    if not cfg.has_option('metadata', 'name'):
        return None
    name = cfg.get('metadata', 'name')
    if cfg.has_option('metadata', 'version'):
        version = cfg.get('metadata', 'version')
    else:
        version = None
    return pkg_resources.Distribution(
               location=path, project_name=name, version=version,
               precedence=pkg_resources.CHECKOUT_DIST)


class LocalSourcesPackageIndex(PackageIndex):
    """
    Like PackageIndex, but can also install packages from local source
    checkouts, the locations of which are added via add_find_links().

    Although PackageIndex supports installing for source distributions on the
    local filesystem, they must be in a tar/zip/etc.  This allows installing
    from an existing source checkout on the local filesystem.
    """

    def process_filename(self, fn, nested=False):
        PackageIndex.process_filename(self, fn, nested)
        dist = distro_from_setup_cfg(fn)
        if dist:
            self.add(dist)

    def fetch_distribution(self, requirement, tmpdir, force_scan=False,
            source=False, develop_ok=False, local_index=None):
        distribute_req = pkg_resources.Requirement.parse('distribute>=0.6.14')
        if pkg_resources.get_distribution('distribute') in distribute_req:
            # The local_index parameter is only in distribute>=0.6.14
            dist = PackageIndex.fetch_distribution(self, requirement, tmpdir,
                                                   force_scan, source,
                                                   develop_ok, local_index)
        else:
            dist = PackageIndex.fetch_distribution(self, requirement, tmpdir,
                                                   force_scan, source,
                                                   develop_ok)
        if dist:
            log.info('Using %s from %s' % (dist, dist.location))
        return dist


class easier_install(easy_install):
    """
    Extension to the easy_install command that uses LocalSourcesPackageIndex as
    its default PackageIndex implementation.
    """

    command_name = 'easy_install'

    create_index = LocalSourcesPackageIndex

    def process_distribution(self, requirement, dist, deps=True, *info):
        """This version of process_distribution will force the package index to
        search for local distributions before going out to PyPI when processing
        a package's dependency.

        It will already do this for the first dependency, but not for
        subsequent dependencies--something I would consider a bug.
        """

        if self.package_index.to_scan is None:
            self.package_index.to_scan = []
        return easy_install.process_distribution(self, requirement, dist, deps,
                                                 *info)
