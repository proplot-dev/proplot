# namedtuple is needed for find_mod_objs so it can have a non-local module

from collections import namedtuple

from ..utils import find_mod_objs


def test_find_mod_objs():
    lnms, fqns, objs = find_mod_objs('sphinx_automodapi')

    # just check for astropy.test ... other things might be added, so we
    # shouldn't check that it's the only thing
    assert lnms == []

    lnms, fqns, objs = find_mod_objs(
        'sphinx_automodapi.tests.test_utils', onlylocals=False)

    assert namedtuple in objs

    lnms, fqns, objs = find_mod_objs(
        'sphinx_automodapi.tests.test_utils', onlylocals=True)
    assert 'namedtuple' not in lnms
    assert 'collections.namedtuple' not in fqns
    assert namedtuple not in objs
