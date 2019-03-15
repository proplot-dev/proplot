# The following tests use a plain Python example module that is at
# sphinx_automodapi.tests.example_module.

# We store different cases in the cases sub-directory of the tests directory

import os
import io
import sys
import glob
import shutil
from itertools import product
from distutils.version import LooseVersion

import pytest

from copy import deepcopy, copy
from sphinx import __version__
from sphinx.util.osutil import ensuredir
from docutils.parsers.rst import directives, roles

SPHINX_LT_16 = LooseVersion(__version__) < LooseVersion('1.6')
SPHINX_LT_17 = LooseVersion(__version__) < LooseVersion('1.7')

if SPHINX_LT_17:
    from sphinx import build_main
else:
    from sphinx.cmd.build import build_main

CASES_ROOT = os.path.join(os.path.dirname(__file__), 'cases')

CASES_DIRS = glob.glob(os.path.join(CASES_ROOT, '*'))

if SPHINX_LT_16 or os.environ.get('TRAVIS_OS_NAME', None) == 'osx':
    PARALLEL = {False}
else:
    PARALLEL = {False, True}


def write_conf(filename, conf):
    with open(filename, 'w') as f:
        for key, value in conf.items():
            f.write("{0} = {1}\n".format(key, repr(conf[key])))


intersphinx_mapping = {
    'python': ('http://docs.python.org/{0}/'.format(sys.version_info[0]), None)
    }

DEFAULT_CONF = {'source_suffix': '.rst',
                'master_doc': 'index',
                'nitpicky': True,
                'extensions': ['sphinx.ext.intersphinx', 'sphinx_automodapi.automodapi'],
                'suppress_warnings': ['app.add_directive', 'app.add_node'],
                'intersphinx_mapping': intersphinx_mapping,
                'nitpick_ignore': [('py:class', 'sphinx_automodapi.tests.example_module.classes.BaseSpam'),
                                   ('py:class', 'sphinx_automodapi.tests.example_module.other_classes.BaseFoo'),
                                   # See the following links for why these classes need to be ignored.
                                   # This only seems to be necessary for Python 2.7.
                                   #
                                   # https://trac.sagemath.org/ticket/19211
                                   # https://stackoverflow.com/q/11417221/3776794
                                   ('py:class', '_abcoll.Sequence'),
                                   ('py:class', '_abcoll.Iterable'),
                                   ('py:class', '_abcoll.Container'),
                                   ('py:class', '_abcoll.Sized')]}


def setup_function(func):
    # This can be replaced with the docutils_namespace context manager once
    # it is in a stable release of Sphinx
    func._directives = copy(directives._directives)
    func._roles = copy(roles._roles)


def teardown_function(func):
    directives._directives = func._directives
    roles._roles = func._roles


@pytest.mark.parametrize(('case_dir', 'parallel'), product(CASES_DIRS, PARALLEL))
def test_run_full_case(tmpdir, case_dir, parallel):

    input_dir = os.path.join(case_dir, 'input')
    output_dir = os.path.join(case_dir, 'output')

    docs_dir = tmpdir.mkdir('docs').strpath

    conf = deepcopy(DEFAULT_CONF)
    conf.update({'automodapi_toctreedirnm': 'api',
                 'automodapi_writereprocessed': True,
                 'automodsumm_writereprocessed': True})

    if os.path.basename(case_dir) in ('mixed_toplevel',
                                      'mixed_toplevel_all_objects'):
        conf['extensions'].append('sphinx_automodapi.smart_resolver')

    start_dir = os.path.abspath('.')

    src_dir = 'src' if 'source_dir' in case_dir else '.'

    ensuredir(os.path.join(docs_dir, src_dir))

    write_conf(os.path.join(os.path.join(docs_dir, src_dir), 'conf.py'), conf)

    for root, dirnames, filenames in os.walk(input_dir):
        for filename in filenames:
            root_dir = os.path.join(docs_dir, os.path.relpath(root, input_dir))
            ensuredir(root_dir)
            input_file = os.path.join(root, filename)
            shutil.copy(input_file, root_dir)

    argv = ['-W', '-b', 'html', src_dir, '_build/html']
    if parallel:
        argv.insert(0, '-j 4')
    if SPHINX_LT_17:
        # As of Sphinx 1.7, the first argument is now no longer ignored
        argv.insert(0, 'sphinx-build')

    try:
        os.chdir(docs_dir)
        status = build_main(argv=argv)
    finally:
        os.chdir(start_dir)

    assert status == 0

    # Check that all expected output files are there and match the reference files
    for root, dirnames, filenames in os.walk(output_dir):
        for filename in filenames:
            path_reference = os.path.join(root, filename)
            path_relative = os.path.relpath(path_reference, output_dir)
            path_actual = os.path.join(docs_dir, path_relative)
            assert os.path.exists(path_actual)
            with io.open(path_actual, encoding='utf8') as f:
                actual = f.read()
            with io.open(path_reference, encoding='utf8') as f:
                reference = f.read()
            assert actual.strip() == reference.strip()
