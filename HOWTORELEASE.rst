Release procedure
=================
ProPlot follows semantic versioning, e.g., v1.0.0. A major version causes incompatible
API changes, a minor version adds functionality, and a patch covers bug fixes.

#. Create a new branch ``release-vX.x.x`` with the version for the release.

  * Update `CHANGELOG.rst` 
  * Make sure all new changes, features are reflected in the documentation.

#. Open a new pull request for this branch targeting `master` 

#. After all tests pass and the PR has been approved, merge the PR into ``master`` 

#. Tag a release and push to github::

    $ git tag -a v1.0.0 -m "Version 1.0.0"
    $ git push origin master --tags

#. Build and publish release on PyPI::

    $ git clean -xfd # remove any files not checked into git
    $ python setup.py sdist bdist_wheel --universal # build package
    $ twine upload dist/* # register and push to pypi
