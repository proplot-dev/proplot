Frequently asked questions
==========================

Why is my text so small?
------------------------
ProPlot makes the default text sizes suitable for publication.


How does ProPlot differ from other matplotlib wrappers?
-------------------------------------------------------

There is already a great matplotlib wrapper called `seaborn <https://seaborn.pydata.org/>`__, and `pandas <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html>`__ and `xarray <http://xarray.pydata.org/en/stable/plotting.html>`__ both offer convenient matplotlib plotting commands. What makes this project different?

While parts of ProPlot were inspired by these projects (in particular, ``rctools.py`` and ``colortools.py`` are modeled after seaborn), the goals are largely different. Seaborn, pandas, and xarray offer tools for generating rigid, simple, nice-looking plots from data stored in special objects (``pandas.DataFrame`` and ``xarray.DataArray``). Similarly, ProPlot uses metadata from these special objects and gives you nice-looking plots out of the box -- but critically, it also permits a *high level of customization*, permits building *complex grids of subplots*, and includes powerful tools for working with *colors* and *geographic datasets*. While seaborn, pandas, and xarray are meant to be used separately from the matplotlib API, ProPlot is built *into the matplotlib API*, thanks to special subclasses of the native matplotlib ``Figure`` and ``Axes`` classes.

In summary, this project is intended to unify the convenience of seaborn, pandas, and xarray plotting with the power and customizability of the underlying matplotlib API.

Why not add to matplotlib directly?
-----------------------------------
Certain parts of ProPlot conflict directly with the matplotlib API. ProPlot enforces a *static* figure layout with the entire subplot grid declared at figure creation time, so that we can implement subplot panels, exert more control on the subplot layout, and replace matplotlib's ``GridSpec`` class with the ``FlexibleGridSpec`` class. By contrast, matplotlib encourages successively adding subplots and panels to existing figures. ProPlot's "smart tight layout" conflicts with matplotlib's `tight layout <https://matplotlib.org/tutorials/intermediate/tight_layout_guide.html>`__ feature by permitting *flexible figure dimensions* to preserve subplot aspect ratios and by permitting *variable subplot spacing* with the ``FlexibleGridSpec`` class.

Other parts of ProPlot would arguably make the matplotlib API redundant if implemented directly. For example, ``Axes.format`` is convenient, but the same tasks can be accomplished with *existing* axes and axis "setter" methods. Also, some of the functionality of ``subplots`` can be replicated with `axes_grid1 <https://matplotlib.org/mpl_toolkits/axes_grid1/index.html>`__. Following `TOOWTDI <https://wiki.python.org/moin/TOOWTDI>`__ philosophy, ProPlot should probably remain here as a separate project.

Nevertheless, if there are any core matplotlib developers reading this, and you think that some of ProPlot's features should be added to matplotlib, please contact me!

