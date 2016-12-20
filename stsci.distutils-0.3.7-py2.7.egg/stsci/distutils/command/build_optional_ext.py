from distutils import log
from distutils.command.build_ext import build_ext
from distutils.errors import DistutilsError, CCompilerError, CompileError
from distutils.util import strtobool

from ConfigParser import ConfigParser


class build_optional_ext(build_ext):
    """This is a version of the build_ext command that allows specific
    extensions to be marked as 'optional'.  If an optional extension fails to
    build, a warning message is displayed but the build does not cancel.

    It should be noted that this functionality already exists in the Python3
    versions of distutils and packaging, but we must provide backwards
    compatibility for it.
    """

    command_name = 'build_ext'

    def _find_optional_extensions(self, extensions):
        """Reads the setup.cfg to determine which extensions were set as
        'optional'.  Optionally, developers may also provide a warning message
        to display (otherwise a default message is used).  This is one way in
        which this command improves on the existing functionality on Python 3.
        """

        # TODO: However, support for the 'optional' extension attribute should
        # be included in d2to1, since it's a legitimate backwards compatibility
        # issue.  But until I do another release it will have to be supported
        # here.

        cfg = ConfigParser()
        try:
            cfg.read('setup.cfg')
        except Exception, e:
            log.warn('Failed to read setup.cfg: %s; proceeding as though '
                     'there are no optional extensions' % e)
            return

        # Map extension names to extension objects
        extensions = dict((ext.name, ext) for ext in extensions)

        for section in cfg.sections():
            if not section.startswith('extension='):
                continue

            # The extension name can be specified with the 'name' option, but
            # if that's missing the name is taken from the section header for
            # now
            if cfg.has_option(section, 'name'):
                name = cfg.get(section, 'name')
            else:
                _, name = section.split('=', 1)

            if name not in extensions:
                # Could happen, for example, if a setup_hook disabled some
                # extensions
                continue

            ext = extensions[name]

            if cfg.has_option(section, 'optional'):
                ext._optional = strtobool(cfg.get(section, 'optional'))
            else:
                ext._optional = False

            if cfg.has_option(section, 'fail_message'):
                ext._fail_message = cfg.get(section, 'fail_message')

    def check_extensions_list(self, extensions):
        build_ext.check_extensions_list(self, extensions)
        self._find_optional_extensions(extensions)

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsError, CompileError), e:
            if not hasattr(ext, '_optional') or not ext._optional:
                raise
            log.warn('building optional extension "%s" failed: %s' %
                     (ext.name, e))
            if hasattr(ext, '_fail_message'):
                log.warn(ext._fail_message)
