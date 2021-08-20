.. _ug_rcmpl: https://matplotlib.org/stable/tutorials/introductory/customizing.html

.. _ug_mplrc: https://matplotlib.org/stable/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files

.. _ug_config:

Configuring ProPlot
===================

Overview
--------

A dictionary-like object named `~proplot.config.rc`, belonging to the
`~proplot.config.Configurator` class, is created when you import ProPlot.
This is your one-stop shop for working with
`matplotlib settings <ug_rcmpl_>`_
stored in `~proplot.config.rc_matplotlib`
(our name for the `~matplotlib.rcParams` dictionary)
and :ref:`ProPlot settings <ug_rcproplot>`
stored in `~proplot.config.rc_proplot`.

To change global settings on-the-fly, simply update `~proplot.config.rc`
using either dot notation or as you would any other dictionary:

.. code-block:: python

  import proplot as pplt
  pplt.rc.name = value
  pplt.rc['name'] = value
  pplt.rc.update(name1=value1, name2=value2)
  pplt.rc.update({'name1': value1, 'name2': value2})

To apply settings to a particular axes or figure, pass the setting
to `proplot.axes.Axes.format` or `proplot.figure.Figure.format`:

.. code-block:: python

  import proplot as pplt
  fig, ax = pplt.subplots()
  ax.format(name1=value1, name2=value2)
  ax.format(rc_kw={'name1': value1, 'name2': value2})

To temporarily modify settings for particular figure(s), pass the setting
to the `~proplot.config.Configurator.context` command:

.. code-block:: python

   import proplot as pplt
   with pplt.rc.context(name1=value1, name2=value2):
       fig, ax = pplt.subplots()
   with pplt.rc.context({'name1': value1, 'name2': value2}):
       fig, ax = pplt.subplots()


In all of these examples, if the setting name contains dots,
you can simply omit the dots. For example, to change the
:rcraw:`title.loc` property, the following approaches are valid:

.. code-block:: python

  import proplot as pplt
  # Apply globally
  pplt.rc.titleloc = value
  pplt.rc.update(titleloc=value)
  # Apply locally
  fig, ax = pplt.subplots()
  ax.format(titleloc=value)

Matplotlib settings
-------------------

Matplotlib settings are natively controlled with the `~matplotlib.rcParams`
dictionary. ProPlot makes this dictionary available in the top-level namespace as
`~proplot.config.rc_matplotlib`. All matplotlib settings can also be changed with
`~proplot.config.rc`. Details on the matplotlib settings can be found on
`this page <ug_rcmpl_>`_.

.. _ug_rcproplot:

ProPlot settings
----------------

.. rubric:: Specific settings

ProPlot-specific settings are controlled with the `~proplot.config.rc_proplot`
dictionary. These settings are not found in `~proplot.config.rc_matplotlib`. They
either control features introduced by ProPlot (e.g., a-b-c labels and geographic
gridlines) or they represent existing matplotlib settings rearranged to be less
verbose or more clearly ordered. Here's a broad overview of the new settings:

* The ``subplots`` category includes settings that control the default
  subplot layout and padding.
* The ``basemap`` setting controls whether basemap is the default geographic plotting
  backend, and the ``cartopy`` category includes cartopy-specific settings.
* The ``abc``, ``title``, and ``label`` categories control a-b-c labels, axes
  titles, and axis labels. The latter two replace ``axes.title`` and ``axes.label``.
* The ``suptitle``, ``leftlabel``, ``toplabel``, ``rightlabel``, and ``bottomlabel``
  categories control the figure titles and subplot row and column labels.
* The ``formatter`` category supersedes matplotlib's ``axes.formatter``
  and includes settings that control the `~proplot.ticker.AutoFormatter` behavior.
* The ``cmap`` category supersedes matplotlib's ``image`` and includes
  settings relevant to colormaps and the `~proplot.colors.DiscreteNorm` normalizer.
* The ``tick`` category supersedes matplotlib's ``xtick`` and ``ytick``
  to simultaneously control *x* and *y* axis tick and tick label settings.
* The matplotlib ``grid`` category includes new settings that control the meridian
  and parallel gridlines and gridline labels managed by `~proplot.axes.GeoAxes`.
* The ``gridminor`` category optionally controls minor gridlines separately
  from major gridlines.
* The ``land``, ``ocean``, ``rivers``, ``lakes``, ``borders``, and ``innerborders``
  categories control geographic content managed by `~proplot.axes.GeoAxes`.

.. rubric:: Meta-settings

Other `~proplot.config.rc` settings may be more accurately described as "meta-settings",
as they change several matplotlib and ProPlot settings at once. For example:

* Setting :rcraw:`font.small` (or, equivalently, :rcraw:`fontsmall`) changes
  the :rcraw:`tick.labelsize`, :rcraw:`grid.labelsize`,
  :rcraw:`legend.fontsize`, and :rcraw:`axes.labelsize`.
* Setting :rcraw:`font.large` (or, equivalently, :rcraw:`texttitlesize`) changes
  the :rcraw:`abc.size`, :rcraw:`title.size`, :rcraw:`suptitle.size`,
  :rcraw:`leftlabel.size`, :rcraw:`toplabel.size`, :rcraw:`rightlabel.size`
  :rcraw:`bottomlabel.size`.
* Setting :rcraw:`meta.color` changes the :rcraw:`axes.edgecolor`,
  :rcraw:`axes.labelcolor` :rcraw:`tick.labelcolor`, :rcraw:`hatch.color`,
  :rcraw:`xtick.color`, and :rcraw:`ytick.color` .
* Setting :rcraw:`meta.width` changes the :rcraw:`axes.linewidth` and the major
  and minor tickline widths :rcraw:`xtick.major.width`, :rcraw:`ytick.major.width`,
  :rcraw:`xtick.minor.width`, and :rcraw:`ytick.minor.width`. The minor tickline widths
  are scaled by :rcraw:`tick.widthratio` (or, equivalently, :rcraw:`tickwidthratio`).
* Setting :rcraw:`tick.len` (or, equivalently, :rcraw:`ticklen`) changes the major and
  minor tickline lengths :rcraw:`xtick.major.size`, :rcraw:`ytick.major.size`,
  :rcraw:`xtick.minor.size`, and :rcraw:`ytick.minor.size`. The minor tickline lengths
  are scaled by :rcraw:`tick.lenratio` (or, equivalently, :rcraw:`ticklenratio`).
* Setting :rcraw:`grid.color`, :rcraw:`grid.linewidth`, :rcraw:`grid.linestyle`,
  or :rcraw:`grid.alpha` also changes the corresponding ``gridminor`` settings. Any
  distinct ``gridminor`` settings must be applied after ``grid`` settings.
* Setting :rcraw:`grid.linewidth` changes the major and minor gridline widths.
  The minor gridline widths are scaled by :rcraw:`grid.widthratio`
  (or, equivalently, :rcraw:`gridwidthratio`).
* Setting :rcraw:`title.border` or :rcraw:`abc.border` to ``True`` automatically
  sets :rcraw:`title.bbox` or :rcraw:`abc.bbox` to ``False``, and vice versa.

.. rubric:: Table of settings

A comprehensive table of the new ProPlot settings is shown below.

.. include:: _static/rctable.rst

.. _ug_proplotrc:

The proplotrc file
------------------

When you import ProPlot for the first time, a ``proplotrc`` file is generated with
all lines commented out. This file is just like `matplotlibrc <ug_mplrc_>`_,
except it controls both matplotlib *and* ProPlot settings. The syntax is essentially
the same as matplotlibrc, and the file path is very similar to matplotlibrc. On most
platforms it is found in ``~/.proplot/proplotrc``, but a loose hidden file in the
home directory named ``~/.proplotrc`` is also allowed (use
`~proplot.config.Configurator.user_file` to print the path). To update this file
after a version change, simply remove it and restart your python session.

To change the global `~proplot.config.rc` settings, edit and uncomment the lines
in the ``proplotrc`` file. To change the settings for a specific project, place a file
named either ``.proplotrc`` or ``proplotrc`` in the same directory as your python
script or jupyter session, or in an arbitrary parent directory. To generate a
``proplotrc`` file containing the settings you have changed during a python session,
use `~proplot.config.Configurator.save` (use `~proplot.config.Configurator.changed`
to preview a dictionary of the changed settings). To explicitly load a ``proplotrc``
file, use `~proplot.config.Configurator.load`.

As an example, a ``proplotrc`` file containing the default settings
is shown below.

.. include:: _static/proplotrc
   :literal:
