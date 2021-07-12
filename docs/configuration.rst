.. _rc_matplotlib: https://matplotlib.org/users/customizing.html

.. _ug_config:

Configuring ProPlot
===================

Overview
--------

A dictionary-like object named `~proplot.config.rc`, belonging to the
`~proplot.config.RcConfigurator` class, is created on import. This is your one-stop
shop for working with `builtin matplotlib global settings <rc_matplotlib_>`_
and the global settings :ref:`added by proplot <rc_proplot>`.
Global settings can be changed on-the-fly using the `~proplot.config.rc`
object:

.. code-block:: python

  import proplot as pplt
  pplt.rc.name = value
  pplt.rc['name'] = value
  pplt.rc.update(name1=value1, name2=value2)
  pplt.rc.update({'name1': value1, 'name2': value2})

To apply settings to a particular axes, pass the setting
to the `~proplot.axes.Axes.format` command:

.. code-block:: python

  import proplot as pplt
  fig, ax = pplt.subplots()
  ax.format(name1=value1, name2=value2)
  ax.format(rc_kw={'name1': value1, 'name2': value2})

To temporarily modify settings for particular figure(s), pass the setting
to the `~proplot.config.RcConfigurator.context` command:

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

Matplotlib settings are natively controlled with the `~matplotlib.rcParams` dictionary.
But all `~matplotlib.rcParams` settings can also be changed with `~proplot.config.rc`.
Details on the matplotlib settings can be found on `this page <rc_matplotlib_>`_.

.. _rc_proplot:

ProPlot settings
----------------

.. rubric:: Specific settings

Some `~proplot.config.rc` settings are not found in `~matplotlib.rcParams`. These
control ProPlot-specific features, like a-b-c labels and geographic content.
Here's a broad overview of the new settings:

* The ``subplots`` category includes settings that control the default
  subplot layout and padding.
* The ``abc``, ``title``, and ``tick`` categories control
  the subplot a-b-c labels, subplot titles, and axis tick labels.
* The ``suptitle``, ``leftlabel``, ``toplabel``, ``rightlabel``,
  and ``bottomlabel`` categories control the figure titles and subplot edge labels.
* The ``formatter`` category is a shorthand for matplotlib's ``axes.formatter``
  and includes new settings that control the default ProPlot axis formatter.
* The matplotlib ``image`` category includes new settings that control the generation
  of automatic colormap levels.
* The matplotlib ``grid`` category includes new settings that control the behavior
  of `~proplot.axes.GeoAxes` gridlines.
* The ``gridminor`` category optionally controls minor gridline settings separately
  from major gridline settings.
* The ``land``, ``ocean``, ``rivers``, ``lakes``, ``borders``, and ``innerborders``
  categories control geographic content managed by `~proplot.axes.GeoAxes`.

.. rubric:: Meta-settings

Other `~proplot.config.rc` settings may be more accurately described as "meta-settings",
as they change several matplotlib and ProPlot settings at once. For example:

* Setting :rcraw:`color` changes the :rcraw:`axes.edgecolor`, :rcraw:`axes.labelcolor`
  :rcraw:`tick.labelcolor`, :rcraw:`hatch.color`, :rcraw:`xtick.color`, and
  :rcraw:`ytick.color` .
* Setting :rcraw:`linewidth` changes the :rcraw:`axes.linewidth`,
  :rcraw:`xtick.major.width`, and :rcraw:`ytick.major.width`.
* Setting :rcraw:`text.labelsize` (or, equivalently, :rcraw:`textlabelsize`) changes
  the :rcraw:`tick.labelsize`, :rcraw:`xtick.labelsize`, :rcraw:`ytick.labelsize`,
  :rcraw:`grid.labelsize`, :rcraw:`legend.fontsize`, and :rcraw:`axes.labelsize`.
* Setting :rcraw:`text.titlesize` (or, equivalently, :rcraw:`texttitlesize`) changes
  the :rcraw:`abc.size`, :rcraw:`title.size`, :rcraw:`suptitle.size`,
  :rcraw:`leftlabel.size`, :rcraw:`toplabel.size`, :rcraw:`rightlabel.size`
  :rcraw:`bottomlabel.size`.
* Setting :rcraw:`tick.len` (or, equivalently, :rcraw:`ticklen`) changes both
  the major tickline and minor tickline length, with the latter scaled
  by :rcraw:`tick.ratio`.
* Setting :rcraw:`grid.color`, :rcraw:`grid.linewidth`, :rcraw:`grid.linestyle`,
  or :rcraw:`grid.alpha` also changes the corresponding ``gridminor`` settings. Any
  distinct ``gridminor`` settings must be applied after ``grid`` settings.

.. rubric:: Table of settings

A comprehensive table of the new ProPlot settings is shown below.

.. include:: _static/rctable.rst

.. _ug_proplotrc:

The proplotrc file
------------------

When you import ProPlot for the first time, a ``.proplotrc`` file is generated and
placed in the home directory. To update this file after a version change, simply
remove it, then import ProPlot again. This file is just like the `matplotlibrc file
<https://matplotlib.org/stable/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files>`__,
but for changing both ProPlot *and* matplotlib settings. The syntax is basically
the same as the ``matplotlibrc`` syntax.

To change the default global settings, simply edit and uncomment the entries listed
in this file. You can also change the settings for individual projects by placing
a file named either ``.proplotrc`` or ``proplotrc`` in the same directory as
your python scripts or jupyter sessions, or in an arbitrary parent directory.
As an example, a ``proplotrc`` file containing the default settings is shown below.

.. include:: _static/proplotrc
   :literal:
