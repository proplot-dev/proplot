.. _ug_config:

Configuring ProPlot
===================

Overview
--------

A special object named `~proplot.config.rc`, belonging to the
`~proplot.config.rc_configurator` class, is created on import. This is your
one-stop shop for changing global settings belonging to any of the following
three categories:

1. `Matplotlib's builtin settings <https://matplotlib.org/users/customizing.html>`__.
   These have the format ``x.y`` or ``x.y.z``.
2. ProPlot's :ref:`added settings <rc_added>`. These also have the format
   ``x.y`` (see below).
3. ProPlot's :ref:`quick settings <rc_quick>`. These have no dots (see below).

You can change settings with the `~proplot.config.rc` object as follows:

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
:rcraw:`title.loc` property, use any of the following approaches:

.. code-block:: python

  import proplot as plot
  # Apply globally
  plot.rc.titleloc = value
  plot.rc.update(titleloc=value)
  # Apply locally
  fig, ax = plot.subplots()
  ax.format(titleloc=value)

.. _rc_params: https://matplotlib.org/users/customizing.html

.. _rc_quick:

Quick settings
--------------

These are *short, simple* names used to change multiple matplotlib and
ProPlot settings at once, as shorthands for settings with longer names, or
for special options. For example, :rcraw:`ticklen` changes the tick length for
the *x* and *y* axes in one go.

================  ==============================================================================================================================================================================================================================================
Key               Description
================  ==============================================================================================================================================================================================================================================
``abc``           Boolean, whether to draw a-b-c labels by default.
``align``         Whether to align axis labels during draw. See `aligning labels <https://matplotlib.org/3.1.1/gallery/subplots_axes_and_figures/align_labels_demo.html>`__.
``alpha``         The opacity of the background axes patch.
``borders``       Boolean, toggles country border lines on and off.
``cmap``          The default colormap.
``coast``         Boolean, toggles coastline lines on and off.
``color``         The color of axis spines, tick marks, tick labels, and labels.
``cycle``         The default color cycle name, used e.g. for lines.
``facecolor``     The color of the background axes patch.
``fontname``      Name of font used for all text in the figure. The default is :ref:`TeX Gyre Heros <ug_fonts>`.
``geogrid``       Boolean, toggles meridian and parallel gridlines on and off.
``grid``          Boolean, toggles major grid lines on and off.
``gridminor``     Boolean, toggles minor grid lines on and off.
``gridratio``     Ratio of minor gridline width to major gridline width.
``inlinefmt``     The inline backend figure format or list thereof. Valid formats include ``'svg'``, ``'pdf'``, ``'retina'``, ``'png'``, and ``jpeg``.
``innerborders``  Boolean, toggles internal border lines on and off, e.g. for states and provinces.
``lakes``         Boolean, toggles lake patches on and off.
``land``          Boolean, toggles land patches on and off.
``large``         Font size for titles, "super" titles, and a-b-c subplot labels.
``linewidth``     Thickness of axes spines and major tick lines.
``lut``           The number of colors to put in the colormap lookup table.
``margin``        The margin of space between axes edges and objects plotted inside the axes, if ``xlim`` and ``ylim`` are unset.
``ocean``         Boolean, toggles ocean patches on and off.
``reso``          Resolution of geographic features, one of ``'lo'``, ``'med'``, or ``'hi'``
``rivers``        Boolean, toggles river lines on and off.
``share``         The axis sharing level, one of ``0``, ``1``, ``2``, or ``3``. See `~proplot.ui.subplots` for details.
``small``         Font size for legend text, tick labels, axis labels, and text generated with `~matplotlib.axes.Axes.text`.
``span``          Boolean, toggles spanning axis labels. See `~proplot.ui.subplots` for details.
``style``         The matplotlib `stylesheet <https://matplotlib.org/3.2.1/gallery/style_sheets/style_sheets_reference.html>`__ name. If ``None``, a custom proplot style is used.
``tickdir``       Major and minor tick direction. Must be one of ``'out'``, ``'in'``, or ``'inout'``.
``ticklen``       Length of major ticks in points.
``ticklenratio``  Ratio of minor tickline length to major tickline length.
``tickpad``       Padding between ticks and tick labels in points.
``titlepad``      Padding between the axes and the title, alias for :rcraw:`axes.titlepad`.
``tickratio``     Ratio of minor tickline width to major tickline width.
``tight``         Boolean, indicates whether to auto-adjust figure bounds and subplot spacings.
================  ==============================================================================================================================================================================================================================================

.. _rc_added:

Added settings
--------------

These are *detailed, descriptive* setting names used to customize things not
covered by `matplotlib's builtin settings
<https://matplotlib.org/users/customizing.html>`__.

The ``subplots`` category controls the default layout for figures and axes.
The ``abc``, ``title``, and ``tick`` categories control a-b-c label, title,
and axis tick label settings. The ``suptitle``, ``leftlabel``, ``toplabel``,
``rightlabel``, and ``bottomlabel`` categories control figure title and edge
label settings.

There are two new additions to the ``image`` category, and the new
``colorbar`` category controls *inset* and *outer*
`~proplot.axes.Axes.colorbar` properties.  The new ``gridminor`` category
controls minor gridline settings, and the new ``geogrid`` category controls
meridian and parallel line settings for `~proplot.axes.GeoAxes`. Note that
when a ``grid`` property is changed, it also changed the corresponding
``gridminor`` property.

Finally, the ``geoaxes``, ``land``, ``ocean``, ``rivers``, ``lakes``,
``borders``, and ``innerborders`` categories control various
`~proplot.axes.GeoAxes` settings. These are used when the boolean toggles for
the corresponding :ref:`quick settings <rc_quick>` are turned on.

===============================  =========================================================================================================================================================================================================================================================
Key                              Description
===============================  =========================================================================================================================================================================================================================================================
``abc.style``                    a-b-c label style. For options, see `~proplot.axes.Axes.format`.
``abc.loc``                      a-b-c label position. For options, see `~proplot.axes.Axes.format`.
``abc.border``                   Boolean, indicates whether to draw a white border around a-b-c labels inside an axes.
``abc.borderwidth``              Width of the white border around a-b-c labels.
``abc.color``                    a-b-c label color.
``abc.size``                     a-b-c label font size.
``abc.weight``                   a-b-c label font weight.
``axes.formatter.zerotrim``      Boolean, indicates whether trailing decimal zeros are trimmed on tick labels.
``axes.formatter.timerotation``  Float, indicates the default *x* axis tick label rotation for datetime tick labels.
``borders.color``                Line color for country borders.
``borders.linewidth``            Line width for country borders.
``bottomlabel.color``            Font color for column labels on the bottom of the figure.
``bottomlabel.size``             Font size for column labels on the bottom of the figure.
``bottomlabel.weight``           Font weight for column labels on the bottom of the figure.
``colorbar.loc``                 Inset colorbar location, options are listed in `~proplot.axes.Axes.colorbar`.
``colorbar.grid``                Boolean, indicates whether to draw borders between each level of the colorbar.
``colorbar.frameon``             Boolean, indicates whether to draw a frame behind inset colorbars.
``colorbar.framealpha``          Opacity for inset colorbar frames.
``colorbar.length``              Length of outer colorbars.
``colorbar.insetlength``         Length of inset colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.width``               Width of outer colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.insetwidth``          Width of inset colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.axespad``             Padding between axes edge and inset colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.extend``              Length of rectangular or triangular "extensions" for panel colorbars. Units are interpreted by `~proplot.utils.units`.
``colorbar.insetextend``         Length of rectangular or triangular "extensions" for inset colorbars. Units are interpreted by `~proplot.utils.units`.
``geoaxes.facecolor``            Face color for the map outline patch.
``geoaxes.edgecolor``            Edge color for the map outline patch.
``geoaxes.linewidth``            Edge width for the map outline patch.
``geogrid.alpha``                Opacity of geographic gridlines.
``geogrid.color``                Color of geographic gridlines.
``geogrid.labelpad``             Default padding in points between map boundary edge and longitude/latitude labels.
``geogrid.labels``               Boolean, indicates whether to label the parallels and meridians.
``geogrid.labelsize``            Font size for latitude and longitude labels. Inherits from ``small``.
``geogrid.latmax``               Absolute latitude in degrees, poleward of which meridian gridlines are cut off.
``geogrid.latstep``              Default interval for parallel gridlines in degrees.
``geogrid.linestyle``            Line style for geographic gridlines.
``geogrid.linewidth``            line width for geographic gridlines.
``geogrid.lonstep``              Default interval for meridian gridlines in degrees.
``geogrid.rotatelabels``         Whether to rotate meridian and parallel gridline labels on cartopy axes.
``gridminor.linewidth``          Minor gridline width.
``gridminor.linestyle``          Minor gridline style.
``gridminor.alpha``              Minor gridline transparency.
``gridminor.color``              Minor gridline color.
``image.levels``                 Default number of levels for ``pcolormesh`` and ``contourf`` plots.
``image.edgefix``                Whether to fix the `white-lines-between-filled-contours <https://stackoverflow.com/q/8263769/4970632>`__ and `white-lines-between-pcolor-rectangles <https://stackoverflow.com/q/27092991/4970632>`__ issues. This slows down figure rendering a bit.
``innerborders.color``           Line color for internal border lines.
``innerborders.linewidth``       Line width for internal border lines.
``land.color``                   Face color for land patches.
``lakes.color``                  Face color for lake patches.
``leftlabel.color``              Font color for row labels on the left-hand side.
``leftlabel.size``               Font size for row labels on the left-hand side.
``leftlabel.weight``             Font weight for row labels on the left-hand side.
``ocean.color``                  Face color for ocean patches.
``rightlabel.color``             Font color for row labels on the right-hand side.
``rightlabel.size``              Font size for row labels on the right-hand side.
``rightlabel.weight``            Font weight for row labels on the right-hand side.
``rivers.color``                 Line color for river lines.
``rivers.linewidth``             Line width for river lines.
``subplots.axwidth``             Default width of each axes. Units are interpreted by `~proplot.utils.units`.
``subplots.panelwidth``          Width of side panels. Units are interpreted by `~proplot.utils.units`.
``subplots.pad``                 Padding around figure edge. Units are interpreted by `~proplot.utils.units`.
``subplots.axpad``               Padding between adjacent subplots. Units are interpreted by `~proplot.utils.units`.
``subplots.panelpad``            Padding between subplots and panels, and between stacked panels. Units are interpreted by `~proplot.utils.units`.
``suptitle.color``               Figure title color.
``suptitle.size``                Figure title font size.
``suptitle.weight``              Figure title font weight.
``tick.labelcolor``              Axis tick label color. Mirrors the *axis* label :rcraw:`axes.labelcolor` setting.
``tick.labelsize``               Axis tick label font size. Mirrors the *axis* label :rcraw:`axes.labelsize` setting.
``tick.labelweight``             Axis tick label font weight. Mirrors the *axis* label :rcraw:`axes.labelweight` setting.
``title.loc``                    Title position. For options, see `~proplot.axes.Axes.format`.
``title.border``                 Boolean, indicates whether to draw a white border around titles inside an axes.
``title.borderwidth``            Width of the white border around titles.
``title.pad``                    The title offset in arbitrary units. Alias for :rcraw:`axes.titlepad`.
``title.color``                  Axes title color.
``title.size``                   Axes title font size.
``title.weight``                 Axes title font weight.
``toplabel.color``               Font color for column labels on the top of the figure.
``toplabel.size``                Font size for column labels on the top of the figure.
``toplabel.weight``              Font weight for column labels on the top of the figure.
===============================  =========================================================================================================================================================================================================================================================

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
