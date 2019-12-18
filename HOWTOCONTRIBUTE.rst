==================
Contribution guide
==================

Contributions are highly welcomed and appreciated.  Every little bit helps,
so please do not hesitate! You can make a high impact on ProPlot just by using it and
reporting `issues <https://github.com/lukelbd/proplot/issues>`__.

The following sections cover some general guidelines
regarding development in ProPlot for maintainers and contributors.
Feel free to suggest improvements or changes in the workflow.

Feature requests and feedback
=============================

We are eager to hear your requests for new features, suggestions regarding the current
API, and so on. You can submit these as
`issues <https://github.com/lukelbd/proplot/issues/new>`__ with the label
"feature."
Please make sure to explain in detail how the feature should work and keep the scope as
narrow as possible. This will make it easier to implement in small pull requests.

If you are feeling inspired, feel free to add the feature yourself!


Report bugs
===========

Bugs should be reported in the `issue tracker <https://github.com/lukelbd/proplot/issues>`__
with the label "bug". When reporting a bug, please include:

* Your operating system name and version, your python version, and your proplot and matplotlib versions.
* If the bug also involves cartopy or basemap, please include these versions as well.
* An example that can be copied and pasted to reproduce the bug.

If you can figure out how to fix the bug, feel free to submit a pull request.

Write tests
===========

Many packages include ``.py`` scripts in a ``tests`` folder
and have the `Travis Continuous Integration <https://travis-ci.com>`__ service
automatically run them. Currently, we do
not use the ``tests`` folder -- we just have Travis run the ``.ipynb`` notebook
examples in the ``docs`` folder (see `.travis.yml`).
However, this is a *major* item on our to-do list!

If you can think of a useful test for ProPlot, feel free to submit a pull request.
Your test will be used in the future.


Write documentation
===================

Documentation can always be improved. For minor changes, you can edit docstrings and documentation files directly in the GitHub web interface without using a local copy.

* The docstrings are written in `reStructuredText <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`__ with `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`__ style headers. They are embedded in the :ref:`API reference` section using a `fork of sphinx-automodapi <https://github.com/lukelbd/sphinx-automodapi>`__. Other sections are written using ``.rst`` and ``.ipynb`` notebook files in the ``docs`` folder. The notebooks are embedded in the User Guide using `nbsphinx <https://nbsphinx.readthedocs.io/en/0.5.0/>`__.
* The `default ReST role <https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-default_role>`__ is ``py:obj``. Please include ``py:obj`` links whenever discussing particular functions or classes -- for example, if you are discussing the `~proplot.axes.Axes.format` method, please write ```~proplot.axes.Axes.format``` rather than ``format``. ProPlot also uses `intersphinx <http://www.sphinx-doc.org/en/stable/ext/intersphinx.html>`__ so you can link to external packages like matplotlib and cartopy.
* When editing the ``.ipynb`` notebook files, make sure to put your example descriptions inside reStructedText cells, not markdown cells. This lets us add sphinx directives and API links to the descriptions. See `this guide <https://nbsphinx.readthedocs.io/en/0.4.3/raw-cells.html#Usage>`__ for how to convert cells to ReST.

To build the documentation locally, use the following commands:

.. code:: bash

   cd docs
   conda env update -f environment.yml
   make html

The built documentation should be available in ``docs/_build/html``.

Preparing pull requests
=======================

#. Fork the
   `proplot GitHub repository <https://github.com/lukelbd/proplot>`__.  It's
   fine to keep "proplot" as the fork repository name because it will live
   under your account.

#. Clone your fork locally using `git <https://git-scm.com/>`__, connect your repository
   to the upstream (main project), and create a branch:

   .. code-block:: bash

      git clone git@github.com:YOUR_GITHUB_USERNAME/proplot.git
      cd proplot
      git remote add upstream git@github.com:lukelbd/proplot.git
      git checkout -b your-branch-name master

   If you need some help with git, follow the
   `quick start guide <https://git.wiki.kernel.org/index.php/QuickStart>`__.

#. Make an editable install of ProPlot by running:

   .. code-block:: bash

      pip install -e .

   This way when you ``import proplot``, the
   local copy is used, rather than the stable version you
   downloaded from PyPi. You can print ``proplot.__file__`` to verify the
   correct version has been imported.

#. Install `pre-commit <https://pre-commit.com>`__ and its hook on the ``proplot`` repo

   .. code-block:: bash

      pip install --user pre-commit
      pre-commit install

   Afterwards ``pre-commit`` will run whenever you commit. https://pre-commit.com/
   is a framework for managing and maintaining multi-language pre-commit hooks to
   ensure code-style and code formatting is consistent.

#. If you intend to make changes or add examples to the ipython notebooks,
   you need to install and configure
   `nbstripout <https://github.com/kynan/nbstripout>`__:

   .. code-block:: bash

      pip install --user nbstripout
      git config --local include.path ../.gitconfig

   This adds the ``proplot/.gitconfig`` file (which is not recognized by git)
   to the local ``proplot/.git/config`` configuration file, which
   defines the filters declared in ``proplot/.gitattributes``. It is necessary
   because git cannot sync repository-specific configuration files.

   After this is done, cell output will be "invisible" to git; the version control
   system only ever "sees" the content written in each cell.
   This makes
   ``git diff``\ s much more legible, significantly reduces the repo size, and
   lets us test notebook examples using
   `nbsphinx <https://nbsphinx.readthedocs.io/en/0.4.3/>`__.

#. You can now edit your local working copy as necessary. Please follow
   the `PEP-8 style guide <https://www.python.org/dev/peps/pep-0008/>`__.
   When committing, ``nbstripout`` will ignore changes to notebook cell output
   and ``pre-commit`` will modify the files as needed, or will generally be clear
   about what you need to do to pass the pre-commit test.

   Please break your edits up into reasonably sized commits:


   .. code-block:: bash

      git commit -a -m "<commit message>"
      git push -u

   The commit messages should be short, sweet, and use the imperative mood,
   e.g. "Fix bug" instead of "Fixed bug".

   ..
      #. Run all the tests. Now running tests is as simple as issuing this command:
         .. code-block:: bash
            coverage run --source proplot -m py.test
         This command will run tests via the ``pytest`` tool against Python 3.7.

#. Create a new changelog entry in ``CHANGELOG.rst``. The entry should be entered as:

   .. code-block::

      <description> (:pr:`<PR number>`) `<author name>`_

   where ``<description>`` is the description of the PR related to the change, ``<PR number>`` is the pull request number, and ``<author name>`` is your first and last name. Add yourself to list of authors at the end of ``CHANGELOG.rst`` if not there, in alphabetical order.

   Make sure to add the changelog entry under one of the valid ``.. rubric:: <heading>`` headings listed at the top of ``CHANGELOG.rst``.

#. Finally, submit a pull request through the GitHub website using this data:

   .. code-block::

      head-fork: YOUR_GITHUB_USERNAME/proplot
      compare: your-branch-name

      base-fork: lukelbd/proplot
      base: master

Note that you can create the pull request while you're working on this. The PR will update
as you add more commits. ProPlot developers and contributors can then review your code
and offer suggestions.


Release procedure
=================

ProPlot follows semantic versioning, e.g. ``vX.Y.Z``. A major version (``X``) causes incompatible
API changes, a minor version (``Y``) adds functionality, and a patch (``Z``) covers bug fixes.

For now, `Luke Davis <https://github.com/lukelbd>`__ is the only one who can publish releases on PyPi, but this will change in the future. Releases should be carried out as follows:


#. Create a new branch ``release-vX.Y.Z`` with the version for the release. In this branch, update ``CHANGELOG.rst``, and make sure all new changes are reflected in the documentation.

   .. code-block:: bash

      git add CHANGELOG.rst
      git commit -m "Changelog updates"


#. Open a new pull request for this branch targeting ``master``.

#. After all tests pass and the pull request has been approved, merge into ``master``.

#. Get the latest version of the master branch:

   .. code-block:: bash

      git checkout master
      git pull

#. Tag the current commit and push to github:

   .. code-block:: bash

      git tag -a vX.Y.Z -m "Version X.Y.Z"
      git push origin master --tags

#. Build and publish release on PyPI:

   .. code-block:: bash

      # Remove previous build products and build the package
      rm -r dist build *.egg-info
      python setup.py sdist bdist_wheel --universal
      # Check the source and upload to the test repository
      twine check dist/*
      twine upload --repository-url https://test.pypi.org/legacy/ dist/*
      # Go to https://test.pypi.org/project/proplot/ and make sure everything looks ok
      # Then make sure the package is installable
      pip install --index-url https://test.pypi.org/simple/ proplot
      # Register and push to pypi
      twine upload dist/*
