Configuration
=============

A special object named ``rc`` (belonging to a class called
``rc_configurator``) is created whenever you import ProPlot. This object
gives you advanced control over the look of your plots. **Use**
``plot.rc`` **as your one-stop shop for changing global settings**.

The rc object
-------------

The ``rc`` object can be used to change built-in
``matplotlib.rcParams`` settings, a few custom “``rcSpecial``” settings,
and some magic “``rcGlobal``” settings that apply to groups of other
settings and keep them synced – e.g., tick, spine, and tick label
colors.

To modify ``rc`` settings, you have three options:

1. Change one global setting using ``plot.rc.name = value`` or ``plot.rc['name'] = value``.
   Note that, for ``rcParams`` settings with ‘dots’ in their name, you will
   have to use ``plot.rc['category.name'] = value``.
2. Update several global settings at once using
   ``plot.rc.update({'name1':value1, 'name2':value2})`` or
   ``plot.rc.update(name1=value1, name2=value2)``, just like you would
   update a dictionary.
3. Change local settings using
   ``ax.format(rc_kw={'name1':value1, 'name2':value2})`` or
   ``ax.format(name1=value1, name2=value2)``. Note that, for this last
   option, **the rc settings will only be applied locally** (i.e. to the
   axes on which ``format()`` is being invoked). This can be convenient for
   (e.g.) highlighting a particular subplot by changing its color.

To **access** a single setting, use ``rc.name`` or ``rc[name]``. To
access a group of setting by category name (e.g., the ``rcParams`` that
look like ``'axes.something'``), use ``rc.axes`` and a **dictionary**
will be returned.

"Global" settings
-----------------

The following is an overview of the available ``rcGlobal`` settings.

* ``color``, ``xcolor``, ``ycolor``: The color of axis spines, tick
  marks, tick labels, and labels. Use the ``x`` and ``y`` versions to just
  change settings for one axis.
* ``facecolor``, ``facehatch``: The
  background color and optional hatching pattern (e.g. ``'xxx'``; see
  `this demo <https://matplotlib.org/gallery/shapes_and_collections/hatch_demo.html>`__)
  for an axes. The latter is useful where you wish to highlight invalid
  (transparent) “NaN” data in a ``pcolormesh`` or ``contourf`` plot.
* ``small``, ``large``: Font size for legend text, tick labels, axis
  labels, and manually placed text created by ``ax.text`` (``small``), and
  for titles, super titles, and a-b-c labels (``large``).
* ``fontname``: The font name used for all text in the figure. ProPlot comes packaged
  with a bunch of desirable fonts, and **changes the default from DejaVu
  Sans (or Bitstream Vera) to Helvetica Neue**. When you first import
  ProPlot, run ``plot.install_fonts()`` and restart your ipython session
  to install them.
* ``linewidth``, ``minorwidth``: thickness of axes spines and major tick
  lines (``linewidth``), and minor tick lines (``minorwidth``, where minor tick line thickness = ``linewidth * minorwidth``).
* ``gridwidth``, ``gridratio``: thickness of gridlines (``gridwidth``), and
  minor gridlines (``gridratio``, where minor gridline thickness
  = ``gridwidth * gridratio``).
* ``gridalpha``, ``gridcolor``, ``gridstyle``: the transparency, color, and line style
  for your major and minor gridlines.
* ``ticklen``, ``tickratio``: length of major ticks (``ticklen``), and
  minor ticks (``tickratio``, where minor tick lengths = ``ticklen * tickratio``).
* ``tickdir``: tick direction, one of ``out``, ``in``, or ``inout``.
* ``abcweight``, ``titleweight``, ``suptitleweight``: the font weight (one of
  ``ultralight``, ``light``, ``normal``, ``medium``, ``demi``, ``bold``,
  ``very bold``, or ``black``) for your title text. Note that many fonts
  only have ``normal`` or ``bold`` available; if you request another
  weight, the “closest” availble weight will be selected.

Note some of these settings can also be controlled using, e.g.,
``ax.format(title_kw={'weight':'bold'})`` instead of
``ax.format(rc_kw={'titleweight':'bold'})``.

Use `~proplot.rcmod.rc_configurator.reset` to reset everything to the initial state. By
default, **settings are reset every time a figure is drawn** -- that is, when
a figure is rendered by the matplotlib backend or saved to file).

