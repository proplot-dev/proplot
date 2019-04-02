# Licensed under a 3-clause BSD style license - see LICENSE.rst

import pytest

from ..utils import iteritems
from . import cython_testpackage

pytest.importorskip('sphinx')  # skips these tests if sphinx not present


class FakeEnv(object):
    """
    Mocks up a sphinx env setting construct for automodapi tests
    """
    def __init__(self, **kwargs):
        for k, v in iteritems(kwargs):
            setattr(self, k, v)


class FakeBuilder(object):
    """
    Mocks up a sphinx builder setting construct for automodapi tests
    """
    def __init__(self, **kwargs):
        self.env = FakeEnv(**kwargs)


class FakeApp(object):
    """
    Mocks up a `sphinx.application.Application` object for automodapi tests
    """
    def __init__(self, srcdir, automodapipresent=True):
        self.builder = FakeBuilder(srcdir=srcdir)
        self.info = []
        self.warnings = []
        self.extensions = []
        if automodapipresent:
            self.extensions.append('sphinx_automodapi.automodapi')

    def info(self, msg, loc):
        self.info.append((msg, loc))

    def warn(self, msg, loc):
        self.warnings.append((msg, loc))


ams_to_asmry_str = """
Before

.. automodsumm:: sphinx_automodapi.automodsumm
    :p:

And After
"""

ams_to_asmry_expected = """\
.. currentmodule:: sphinx_automodapi.automodsumm

.. autosummary::
    :p:

    Automoddiagram
    Automodsumm
    automodsumm_to_autosummary_lines
    generate_automodsumm_docs
    process_automodsumm_generation
"""


def test_ams_to_asmry(tmpdir):
    from ..automodsumm import automodsumm_to_autosummary_lines

    fi = tmpdir.join('automodsumm.rst')
    fi.write(ams_to_asmry_str)

    fakeapp = FakeApp(srcdir='')
    resultlines = automodsumm_to_autosummary_lines(str(fi), fakeapp)

    assert '\n'.join(resultlines) == ams_to_asmry_expected


ams_cython_str = """
Before

.. automodsumm:: apyhtest_eva.unit02
    :functions-only:
    :p:

And After
"""

ams_cython_expected = """\
.. currentmodule:: apyhtest_eva.unit02

.. autosummary::
    :p:

    pilot
"""


def test_ams_cython(tmpdir, cython_testpackage):
    from ..automodsumm import automodsumm_to_autosummary_lines

    fi = tmpdir.join('automodsumm.rst')
    fi.write(ams_cython_str)

    fakeapp = FakeApp(srcdir='')
    resultlines = automodsumm_to_autosummary_lines(str(fi), fakeapp)

    assert '\n'.join(resultlines) == ams_cython_expected
