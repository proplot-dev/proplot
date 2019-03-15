import os
import sys
import subprocess as sp

from textwrap import dedent

import pytest


@pytest.fixture
def cython_testpackage(tmpdir, request):
    """
    Creates a trivial Cython package for use with tests.
    """

    test_pkg = tmpdir.mkdir('test_pkg')
    test_pkg.mkdir('apyhtest_eva').ensure('__init__.py')
    test_pkg.join('apyhtest_eva').join('unit02.pyx').write(dedent("""\
        def pilot():
            \"\"\"Returns the pilot of Eva Unit-02.\"\"\"

            return True
    """))

    import sphinx_automodapi  # noqa

    test_pkg.join('setup.py').write(dedent("""\
        import sys

        sys.path.insert(0, {0!r})

        from os.path import join
        from setuptools import setup, Extension

        NAME = 'apyhtest_eva'
        VERSION = 0.1
        RELEASE = True


        setup(
            name=NAME,
            version=VERSION,
            ext_modules=[Extension('apyhtest_eva.unit02',
                                   [join('apyhtest_eva', 'unit02.pyx')])]
        )
    """.format(os.path.dirname(sphinx_automodapi.__path__[0]))))

    test_pkg.chdir()
    # Build the Cython module in a subprocess; otherwise strange things can
    # happen with Cython's global module state
    sp.call([sys.executable, 'setup.py', 'build_ext', '--inplace'])

    sys.path.insert(0, str(test_pkg))
    import apyhtest_eva.unit02

    def cleanup(test_pkg=test_pkg):
        for modname in ['apyhtest_eva', 'apyhtest_eva.unit02']:
            try:
                del sys.modules[modname]
            except KeyError:
                pass

        sys.path.remove(str(test_pkg))

    request.addfinalizer(cleanup)

    return test_pkg
