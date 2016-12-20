"""AST Visitors

Currently only uses one for collecting a list of import statements.
Unfortunately two versions of this have to be implemented: One for 2.6 and up
and a different version for 2.5.
"""


import os
import sys

from distutils import log

try:
    import ast # Python >= 2.6

    def walk(filename, visitor):
        """Generate an AST for the given filename and walk over it using
        the given visitor instance.
        """

        filename = os.path.abspath(filename)

        try:
            tree = ast.parse(open(filename, 'r').read())
        except SyntaxError:
            if sys.version_info[0] < 3:
                e = sys.exc_info()[1]
                log.warn('SyntaxError while parsing file %s: %s' %
                         (filename, str(e)))
                return
            # We're probably in Python 3 and looking at a file intended for
            # Python 2.  Otherwise there's an unintended SyntaxError in the
            # file, so there are bigger problems anyways
            try:
                import lib2to3.refactor

                tool = StringRefactoringTool(
                    lib2to3.refactor.get_fixers_from_package('lib2to3.fixes'))
                tool.refactor_file(filename, write=True)
                tree = ast.parse(tool.refactored[filename])
            except ImportError:
                # Without 2to3 we can't do much more.
                # TODO: Issue a warning?
                return

        visitor.visit(tree)


    class ImportVisitor(ast.NodeVisitor):
        def __init__(self):
            self.imports = set()
            self.importfroms = set()

        def visit_Import(self, node):
            for name in node.names:
                self.imports.add((name.name, name.asname))

        def visit_ImportFrom(self, node):
            for name in node.names:
                self.importfroms.add((node.module, name.name, name.asname))

except ImportError:
    import compiler

    def walk(filename, visitor):
        tree = compiler.parseFile(filename)
        compiler.walk(tree, visitor)

    class ImportVisitor(compiler.visitor.ASTVisitor):
        def __init__(self):
            self.imports = set()
            self.importfroms = set()

        def visitImport(self, node):
            for name in node.names:
                self.imports.add(name)

        def visitFrom(self, node):
            for name in node.names:
                self.importfroms.add((node.modname, name[0], name[1]))


if sys.version_info[0] >= 3:
    try:
        import lib2to3.refactor

        class StringRefactoringTool(lib2to3.refactor.RefactoringTool):
            """A RefactoringTool that saves refactored files as strings in the
            self.refactored dict rather than outputting to actual files.

            This is used in case we're running in Python 3 and need to refactor
            a file before parsing its syntax tree.
            """

            def __init__(self, fixer_names, options=None, explicit=None):
                super(StringRefactoringTool, self).__init__(fixer_names,
                                                            options,
                                                            explicit)
                self.refactored = {}

            def write_file(self, new_text, filename, old_text, encoding=None):
                self.refactored[filename] = new_text

    except ImportError:
        pass
