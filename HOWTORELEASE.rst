=================
Release procedure
=================

ProPlot follows semantic versioning, e.g. v1.0.0. A major version causes incompatible
API changes, a minor version adds functionality, and a patch covers bug fixes.

#. Create a new branch ``release-vX.x.x`` with the version for the release.

#. Update ``CHANGELOG.rst``, and make sure all new changes, features are reflected in the documentation.

#. Open a new pull request for this branch targeting ``master``.

#. After all tests pass and the PR has been approved, merge the PR into ``master``.

#. Pull down the new version of master:

   .. code-block:: bash

      git checkout master
      git pull

#. Tag a release and push to github:

   .. code-block:: bash

      git tag -a v1.0.0 -m "Version 1.0.0"
      git push origin master --tags

#. Build and publish release on PyPI:

   .. code-block:: bash

      rm -r dist build *.egg-info # remove previous build products
      python setup.py sdist bdist_wheel --universal # build package
      twine check dist/* # check that the README is valid
      twine upload --repository-url https://test.pypi.org/legacy/ dist/* # test
      pip install --index-url https://test.pypi.org/simple/ proplot # test
      twine upload dist/* # register and push to pypi

