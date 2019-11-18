==================
Contribution Guide
==================

Contributions are highly welcomed and appreciated.  Every little help counts,
so do not hesitate! You can make a high impact on ``proplot`` just by using it and
reporting `issues <https://github.com/lukelbd/proplot/issues>`_.

The following sections cover some general guidelines
regarding development in ``proplot`` for maintainers and contributors.
Nothing here is set in stone and can't be changed.
Feel free to suggest improvements or changes in the workflow.

Thanks to Riley Brady and his `climpred <https://github.com/lukelbd/proplot/blob/master/HOWTOCONTRIBUTE.rst>`_ project for providing a template for this guide.

.. contents:: Contribution links
   :depth: 2


Feature requests and feedback
=============================

We are eager to hear about your requests for new features and any suggestions about the
API, infrastructure, and so on. Feel free to submit these as
`issues <https://github.com/lukelbd/proplot/issues/new>`_ with the label "feature request."

Please make sure to explain in detail how the feature should work and keep the scope as
narrow as possible. This will make it easier to implement in small PRs.


.. _reportbugs:

Report bugs
===========

Report bugs for ``proplot`` in the `issue tracker <https://github.com/lukelbd/proplot/issues>`_
with the label "bug".

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting,
  specifically the Python interpreter version, installed libraries, and ``proplot``
  version.
* Detailed steps to reproduce the bug.

If you can write a demonstration test that currently fails but *should* pass,
that is a very useful commit to make as well, even if you cannot fix the bug itself.


.. _fixbugs:

Fix bugs
========

Look through the `GitHub issues for bugs <https://github.com/lukelbd/proplot/labels/bug>`_.

Talk to developers to find out how you can fix specific bugs.

Write documentation
===================

``proplot`` could always use more documentation.  What exactly is needed?

* More complementary documentation.  Have you perhaps found something unclear?
* Docstrings.  There can never be too many of them.
* Example notebooks with different Earth System Models, lead times, etc. -- they're all very appreciated.

You can also edit documentation files directly in the GitHub web interface,
without using a local copy.  This can be convenient for small fixes.

Our documentation is written in reStructuredText. You can follow our conventions in already written documents. Some helpful guides are located `here <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_ and `here <https://github.com/ralsina/rst-cheatsheet/blob/master/rst-cheatsheet.rst>`_.

.. note::
    Build the documentation locally with the following command:

    .. code:: bash

        $ conda env update -f ci/environment-dev-3.6.yml
        $ cd docs
        $ make html

    The built documentation should be available in the ``docs/build/``.

If you need to add new functions to the API, run ``sphinx-autogen -o api api.rst`` from the ``docs/source`` directory and add the functions to ``api.rst``.

 .. _`pull requests`:
.. _pull-requests:


Preparing Pull Requests
=======================

.. _submitfeedback:

#. If you intend to make changes / add examples to the ipython notebooks,
   you need to install [nbstripout](https://github.com/kynan/nbstripout)
   with ``pip install nbstripout``. This deletes notebook cell output so
   `nbsphinx <https://nbsphinx.readthedocs.io/en/0.4.3/>`__ always reruns them
   after git pushes.

#. Fork the
   `proplot GitHub repository <https://github.com/lukelbd/proplot>`__.  It's
   fine to use ``proplot`` as your fork repository name because it will live
   under your user.

#. Clone your fork locally using `git <https://git-scm.com/>`_, connect your repository
   to the upstream (main project), and create a branch::

    $ git clone git@github.com:YOUR_GITHUB_USERNAME/proplot.git
    $ cd proplot
    $ git remote add upstream git@github.com:lukelbd/proplot.git

    # now, to fix a bug or add feature create your own branch off "master":

    $ git checkout -b your-bugfix-feature-branch-name master

   If you need some help with Git, follow this quick start
   guide: https://git.wiki.kernel.org/index.php/QuickStart

#. Install dependencies into a new conda environment::

    $ conda env update -f ci/environment-dev-3.7.yml
    $ conda activate proplot-dev

#. Make an editable install of proplot by running::

    $ pip install -e .

#. Break your edits up into reasonably sized commits.

    $ git commit -a -m "<commit message>"
    $ git push -u

#. Run all the tests

   Now running tests is as simple as issuing this command::

    $ coverage run --source proplot -m py.test


   This command will run tests via the "pytest" tool against Python 3.6.


#. Create a new changelog entry in ``CHANGELOG.rst``:

   - The entry should be entered as:

    <description> (``:pr:`#<pull request number>```) ```<author's names>`_``

    where ``<description>`` is the description of the PR related to the change and ``<pull request number>`` is
    the pull request number and ``<author's names>`` are your first and last names.

   - Add yourself to list of authors at the end of ``CHANGELOG.rst`` file if not there yet, in alphabetical order.

 #. Add yourself to the `contributors <https://proplot.readthedocs.io/en/latest/contributors.html>_` list via ``docs/source/contributors.rst``.

#. Finally, submit a pull request through the GitHub website using this data::

    head-fork: YOUR_GITHUB_USERNAME/proplot
    compare: your-branch-name

    base-fork: lukelbd/proplot
    base: master

Note that you can create the Pull Request while you're working on this. The PR will update
as you add more commits. ``proplot`` developers and contributors can then review your code
and offer suggestions.
