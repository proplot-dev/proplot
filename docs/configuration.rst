.. _rc_matplotlib: https://matplotlib.org/users/customizing.html

.. _ug_config:

Configuring ProPlot
===================

Overview
--------

A special object named `~proplot.config.rc`, belonging to the
`~proplot.config.RcConfigurator` class, is created on import. This is your one-stop
shop for working with `builtin matplotlib global settings <rc_matplotlib>`_
and the global settings :ref:`added by proplot <rc_proplot>`.
Global settings can be changed on-the-fly using the `~proplot.config.rc`
object as follows:

.. code-block:: python

  import proplot as plot
  plot.rc.name = value
  plot.rc['name'] = value
  plot.rc.update(name1=value1, name2=value2)
  plot.rc.update({'name1': value1, 'name2': value2})

To apply settings to a particular axes, pass the setting
to the `~proplot.axes.Axes.format` command using either
of the following approaches:

.. code-block:: python

  import proplot as plot
  fig, ax = plot.subplots()
  ax.format(name1=value1, name2=value2)
  ax.format(rc_kw={'name1': value1, 'name2': value2})

In all of these examples, if the setting name contains dots,
you can simply omit the dots. For example, to change the
:rcraw:`title.loc` property, the following approaches are valid:

.. code-block:: python

  import proplot as plot
  # Apply globally
  plot.rc.titleloc = value
  plot.rc.update(titleloc=value)
  # Apply locally
  fig, ax = plot.subplots()
  ax.format(titleloc=value)

Matplotlib settings
-------------------

Details on the matplotlib settings can be found on `this page <rc_matplotlib_>`_.

.. _rc_proplot:

ProPlot settings
----------------

ProPlot adds several settings to customize things not covered by
`matplotlib's builtin settings <https://matplotlib.org/users/customizing.html>`__.
Some of these settings may be more accurately described as "meta-settings",
as they change several matplotlib settings at once (for example,
:rcraw:`linewidth` changes axes edge widths, gridline widths, and tick widths).
Other settings are for specific features controlled by ProPlot, like a-b-c labels.

The ``subplots`` category controls the default layout for figures and axes.
The ``abc``, ``title``, and ``tick`` categories control a-b-c label, title,
and axis tick label settings. The ``suptitle``, ``leftlabel``, ``toplabel``,
``rightlabel``, and ``bottomlabel`` categories control figure title and edge
label settings. There are two new additions to the ``image`` category, and the
new ``colorbar`` category controls *inset* and *outer*
`~proplot.axes.Axes.colorbar` properties. There is also a new ``gridminor``
category for minor gridline settings (note that ``gridminor`` inherits
``grid`` properties when they are changed).
Finally, the ``land``, ``ocean``, ``rivers``, ``lakes``,
``borders``, and ``innerborders`` categories control various
`~proplot.axes.GeoAxes` settings.

.. include:: _static/rctable.rst

.. _ug_proplotrc:

The .proplotrc file
-------------------

When you install ProPlot for the first time, a ``.proplotrc`` file is generated
and placed in your home directory. This is just like the `matplotlibrc file\
<https://matplotlib.org/3.1.1/tutorials/introductory/customizing.html#customizing-with-matplotlibrc-files>`__,
but for changing both ProPlot *and* matplotlib settings. The syntax is basically
the same as the ``matplotlibrc`` syntax.

To change the default global settings, simply edit and uncomment the entries
listed in this file. You can also change the settings for individual projects
by placing a ``.proplotrc`` file in the same directory as your python scripts
or jupyter notebooks, or in an arbitrary parent directory. As an example,
a ``.proplotrc`` file containing the default settings is shown below.

.. include:: _static/proplotrc
   :literal:
