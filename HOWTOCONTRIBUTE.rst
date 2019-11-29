==================
Contribution guide
==================

Contributions are highly welcomed and appreciated.  Every little help counts,
so do not hesitate! You can make a high impact on ProPlot just by using it and
reporting `issues <https://github.com/lukelbd/proplot/issues>`__.

The following sections cover some general guidelines
regarding development in ProPlot for maintainers and contributors.
Nothing here is set in stone and can't be changed.
Feel free to suggest improvements or changes in the workflow.

Feature requests and feedback
=============================

We are eager to hear about your requests for new features and any suggestions about the
API, infrastructure, and so on. Feel free to submit these as
`issues <https://github.com/lukelbd/proplot/issues/new>`__ with the label "feature."

Please make sure to explain in detail how the feature should work and keep the scope as
narrow as possible. This will make it easier to implement in small PRs.


Report bugs
===========

Report bugs for ProPlot in the `issue tracker <https://github.com/lukelbd/proplot/issues>`__
with the label "bug".

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting,
  specifically the Python interpreter version, installed libraries, and ProPlot
  version.
* Detailed steps to reproduce the bug.

If you can write a demonstration test that currently fails but *should* pass,
that is a very useful commit to make as well, even if you cannot fix the bug itself.


Fix bugs
========

Look through the `GitHub issues for bugs <https://github.com/lukelbd/proplot/labels/bug>`__.

Talk to developers to find out how you can fix specific bugs.

Write documentation
===================

ProPlot could always use better documentation. For small changes, you can edit documentation files directly in the GitHub web interface,
without using a local copy.

* The documentation is written in reStructuredText with `numpydoc <https://numpydoc.readthedocs.io/en/latest/>`__ style headers.
* The `default ReST role <https://www.sphinx-doc.org/en/master/usage/configuration.html#confval-default_role>`__ is ``'py:obj'``. This is meant to encourage populating docstrings with links to the API reference. ProPlot uses `intersphinx <http://www.sphinx-doc.org/en/stable/ext/intersphinx.html>`__, so you can also link to sphinx documentation from other projects, e.g. matplotlib. In general, you should compress your links with a tilde, e.g. ``~path.to.function``.
* When editing the ipython notebooks found in ``docs``, make sure to put your example descriptions inside reStructedText cells, not markdown cells. This lets us populate the descriptions with sphinx links. See `this guide <https://nbsphinx.readthedocs.io/en/0.4.3/raw-cells.html#Usage>`__ for how to convert cells to ReST.

Some helpful ReST guides are located `here <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`__ and `here <https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst>`__.

.. note::
    To build the documentation locally, use the following commands:

    .. code:: bash

        cd docs
        pip install requirements.txt
        make html

    The built documentation should be available in the ``docs/_build/html``.

Preparing pull requests
=======================

#. Fork the
   `proplot GitHub repository <https://github.com/lukelbd/proplot>`__.  It's
   fine to use ProPlot as your fork repository name because it will live
   under your user.

#. Clone your fork locally using `git <https://git-scm.com/>`__, connect your repository
   to the upstream (main project), and create a branch:

   .. code-block:: bash

      git clone git@github.com:YOUR_GITHUB_USERNAME/proplot.git
      cd proplot
      git remote add upstream git@github.com:lukelbd/proplot.git
      git checkout -b your-bugfix-feature-branch-name master

   If you need some help with git, follow the
   `quick start guide <https://git.wiki.kernel.org/index.php/QuickStart>`__.

#. Install `pre-commit <https://pre-commit.com>`_ and its hook on the ``proplot`` repo

   .. code-block:: bash

      pip install --user pre-commit
      pre-commit install

   Afterwards ``pre-commit`` will run whenever you commit. https://pre-commit.com/
   is a framework for managing and maintaining multi-language pre-commit hooks to
   ensure code-style and code formatting is consistent.

   You can now edit your local working copy as necessary. Please follow
   PEP-8 naming conventions. When committing, ``pre-commit`` will modify the
   files as needed, or will generally be clear about what you need to do to
   pass the commit test.

#. If you intend to make changes / add examples to the ipython notebooks,
   you need to install and configure
   `nbstripout <https://github.com/kynan/nbstripout>`__ with

   .. code-block:: bash

      pip install --user nbstripout
      git config --local include.path ../.gitconfig

   This strips notebook cell output when files are staged, which reduces the
   repo storage size and lets us use
   `nbsphinx <https://nbsphinx.readthedocs.io/en/0.4.3/>`__
   to test each ``git push``.

   The ``git config`` command associates the filters declared in
   ``proplot/.gitattributes`` with the operations described in ``proplot/.gitconfig``
   by adding them to the *recognized* local configuration file
   ``proplot/.git/config``.

#. Make an editable install of ProPlot by running:

   .. code-block:: bash

      pip install -e .

   This way when you ``import proplot``, your
   local copy is used. You can print ``proplot.__file__`` to verify this.
   Make sure matplotlib is already installed.

#. Break your edits up into reasonably sized commits.

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

#. Create a new changelog entry in ``CHANGELOG.rst``:

   - The entry should be entered as:

     .. code-block::

         <description> (:pr:`<PR number>`) `<author name>`_

     where ``<description>`` is the description of the PR related to the change, ``<PR number>`` is the pull request number, and ``<author name>`` is your first and last name.

   - Add yourself to list of authors at the end of ``CHANGELOG.rst`` file if not there yet, in alphabetical order.

#. Finally, submit a pull request through the GitHub website using this data:

   .. code-block::

      head-fork: YOUR_GITHUB_USERNAME/proplot
      compare: your-branch-name

      base-fork: lukelbd/proplot
      base: master

Note that you can create the Pull Request while you're working on this. The PR will update
as you add more commits. ProPlot developers and contributors can then review your code
and offer suggestions.

