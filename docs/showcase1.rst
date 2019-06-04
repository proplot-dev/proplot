Introduction
============

Matplotlib APIs
---------------

Matplotlib has two APIs – the “pyplot” API (which is MATLAB-like), and
the “object-oriented” API (which is more “pythonic”, more clear, more
flexible, and you should consider using!).

Contrary to the similar names, this package is not meant to be a pyplot
replacement – it adds to the “object-oriented” API by subclassing
matplotlib Artists classes, like `~matplotlib.axes.Axes` and
`~matplotlib.figure.Figure`.

Since `~matplotlib.pyplot` mostly just calls the object-oriented API
under-the-hood, you can still use it to access some of ProPlot’s
features – but this is discouraged.

.. code:: ipython3

    import matplotlib.pyplot as plt
    import numpy as np
    plt.figure(figsize=(5,3))
    plt.plot(np.random.rand(10,10), lw=3)
    plt.title('PyPlot API (discouraged)')
    plt.xlabel('x axis')
    plt.ylabel('y axis')







.. image:: showcase/showcase_2_1.svg


.. code:: ipython3

    import matplotlib.pyplot as plt
    import numpy as np
    f, ax = plt.subplots(figsize=(5,3))
    ax.plot(np.random.rand(10,10), lw=3)
    ax.set_title('Object-oriented API (recommended)')
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')







.. image:: showcase/showcase_3_1.svg


The subplots command
--------------------

The `~proplot.subplots.subplots` command is your gateway to all of
ProPlot’s features. Its usage is sort of like the pyplot
`~matplotlib.pyplot.subplots` version, but it is packed with new
features and generates a subclassed figure and specially subclassed
axes.

.. code:: ipython3

    import proplot as plot
    f, ax = plot.subplots(width=2)
    ax.format(title='ProPlot API', xlabel='x axis', ylabel='y axis')
    f, axs = plot.subplots(ncols=3, nrows=2, width=5)
    axs.format(title='Axes', suptitle='ProPlot API', xlabel='x axis', ylabel='y axis')



.. image:: showcase/showcase_6_0.svg



.. image:: showcase/showcase_6_1.svg


Most matplotlib sizing arguments assume the units “inches” or some
“relative” unit size – e.g. relative to the axes width. With ProPlot,
most sizing arguments are interpreted as follows: if numeric, the units
are inches, and if string, the units are interpreted by
`~proplot.utils.units` (see `~proplot.utils.units` documentation for
a table). Note this means even `~matplotlib.gridspec.GridSpec`
arguments like ``wspace`` and ``hspace`` accept physical units (see
`~proplot.subplots.subplots` for details). I recognize that the rest
of the world doesn’t use “inches”, so I thought this would be useful.

.. code:: ipython3

    import proplot as plot
    f, axs1 = plot.subplots(ncols=2, axwidth=1, height='45mm')
    f, axs2 = plot.subplots(width='5cm', aspect=(2,1))
    f, axs3 = plot.subplots(height='150pt', aspect=0.5)
    for axs in (axs1,axs2,axs3):
        axs.format(suptitle='Title', xlabel='x axis', ylabel='y axis')



.. image:: showcase/showcase_8_0.svg



.. image:: showcase/showcase_8_1.svg



.. image:: showcase/showcase_8_2.svg


Subplot labelling is another useful ProPlot feature. The label order is
row-major by default; to change this, use the
`~proplot.subplots.subplots` ``order`` keyword arg. Change the label
position with the ``abcpos`` `~proplot.rcmod` option, or the label
style with the ``abcformat`` `~proplot.rcmod` option. Toggle labelling
with ``abc=True``. See :ref:`The format command` and
:ref:`Rc settings control` for details.

.. code:: ipython3

    import proplot as plot
    f, axs = plot.subplots(nrows=2, ncols=2, order='F', axwidth=1.5)
    axs.format(abc=True, abcpos='ol', abcformat='A.', xlabel='x axis', ylabel='y axis', suptitle='Subplots with column-major labelling')
    f, axs = plot.subplots(nrows=8, ncols=8, axwidth=0.5, flush=True) # not 
    axs.format(abc=True, abcpos='ir', xlabel='x axis', ylabel='y axis', xticks=[], yticks=[], suptitle='Grid of "flush" subplots')



.. image:: showcase/showcase_10_0.png
   :width: 364px
   :height: 393px



.. image:: showcase/showcase_10_1.png
   :width: 562px
   :height: 572px


.. code:: ipython3

    import proplot as plot
    f, axs = plot.subplots(nrows=2, ncols=2, order='F', axwidth=1.5)
    axs.format(abc=True, abcpos='ol', abcformat='A.', xlabel='x axis', ylabel='y axis', suptitle='Subplots with column-major labelling')
    f, axs = plot.subplots(nrows=8, ncols=8, axwidth=0.5, flush=True) # not 
    axs.format(abc=True, abcpos='ir', xlabel='x axis', ylabel='y axis', xticks=[], yticks=[], suptitle='Grid of "flush" subplots')



.. image:: showcase/showcase_11_0.svg



.. image:: showcase/showcase_11_1.svg


To set up a complex grid of subplots, use a 2D array of integers. You
can think of this array as a “picture” of your figure. This lets you
build the below grid in just one line of code, instead of 6 lines. The
numbering determines the order of a-b-c labels. See
`~proplot.subplots.subplots` for details.

.. code:: ipython3

    # Arbitrarily complex array of subplots, with shared/spanning x/y axes detected automatically
    import proplot as plot
    import numpy as np
    f, axs = plot.subplots([[1, 1, 2], [1, 1, 6], [3, 4, 4], [3, 5, 5]], span=1, share=3, width=5)
    axs.format(suptitle='Complex subplot grid with axis-sharing + spanning labels', xlabel='time (seconds)', ylabel='temperature (K)', abc=True)
    axs[0].plot(2*(np.random.rand(100,5)-0.5).cumsum(axis=0), lw=2)







.. image:: showcase/showcase_13_1.svg


Automatic subplot spacing
-------------------------

With ProPlot, you will always get just the right amount of spacing
between subplots so that elements don’t overlap, and just the right
amount of space around the figure edge so that labels and whatnot are
not cut off. Furthermore, despite all of the complex adjustments this
requires, the original subplot aspect ratios are **always preserved**.
Even when axes panels are present, the main subplot aspect ratios will
stay fixed (see below for more on panels).

You can disable this feature by passing ``tight=False`` to
`~proplot.subplots.subplots`, but it is unbelievably useful. It works
by scaling either the figure width or height dimension (whichever one
you didn’t specify) such that the subplot aspect ratios will not change,
and by taking advantage of ProPlot’s subplot layout restrictions. Some
examples are below.

Sometimes, ``tight=True`` is not possible (when using the cartopy
``set_extent`` method or when using cartopy meridian and parallel
labelling; a warning will be raised in these instances). Even when
``tight=False``, ProPlot tries to make the default spacing reasonable.

.. code:: ipython3

    import proplot as plot
    for share in (3,0):
        f, axs = plot.subplots(nrows=3, ncols=3, aspect=1, axwidth=1, share=share, span=False, tight=True)
        axs[4].format(title='title\ntitle\ntitle', suptitle='"Smart tight layout" automatic spacing')
        axs[1].format(ylabel='ylabel', xlabel='xlabel')



.. image:: showcase/showcase_16_0.svg



.. image:: showcase/showcase_16_1.svg


.. code:: ipython3

    import proplot as plot
    f, axs = plot.subplots([[1,2],[3,2],[3,4]], share=0, span=0, axwidth=1.5)
    axs[0].format(xlabel='xlabel\nxlabel\nxlabel', title='Title', suptitle='"Smart tight layout" automatic spacing')
    axs[1].format(ylabel='ylabel\nylabel', xformatter='null', yticklabelloc='both')
    axs[2].format(yformatter='null', title='Title', ytickloc='both')
    axs[3].format(yformatter='null', xlabel='xlabel\nxlabel\nxlabel')



.. image:: showcase/showcase_17_0.svg


.. code:: ipython3

    import proplot as plot
    f, axs = plot.subplots(axwidth=2, ncols=2, span=False, share=0, axpanels='lr', axpanels_kw={'rshare':False})
    axs.format(ylabel='ylabel', xlabel='xlabel')
    axs[0].lpanel.format(ytickloc='right', yticklabelloc='right')
    axs[0].rpanel.format(ylabel='ylabel', ytickloc='right', yticklabelloc='right',
                         suptitle='"Smart tight layout" automatic spacing with panels', collabels=['Column 1', 'Column 2'])



.. image:: showcase/showcase_18_0.svg


The format command
------------------

The `~proplot.subplots.subplots` method populates the
`~proplot.subplots.Figure` object with either `~proplot.axes.XYAxes`
(for cartesian axes) or `~proplot.axes.MapAxes` (for cartopy or
basemap map projection axes). Both of these classes inherit from the
base class `~proplot.axes.BaseAxes`.

The **most important** new method you need to know is
`~proplot.axes.BaseAxes.format`. This is your one-stop-shop for
changing axis labels, tick labels, titles, etc. Keyword args passed to
this function are interpreted as follows:

1. Any keyword arg matching the name of a ProPlot or native matplotlib
   “rc” setting will be applied to the axes. If the name has “dots”,
   simply omit them. See the `~proplot.rcmod` documentation for
   details.
2. Remaining keyword args are passed to the ``smart_update`` methods of
   the top-level class – that is, the `~proplot.axes.XYAxes`
   `~proplot.axes.XYAxes.smart_update` or `~proplot.axes.MapAxes`
   `~proplot.axes.MapAxes.smart_update` methods. Use these to change
   settings specific to Cartesian axes or specific to map projections,
   like tick locations and toggling geographic features.
3. Finally, the remaining keyword args are passed to the
   `~proplot.axes.BaseAxes` `~proplot.axes.BaseAxes.smart_update`
   method. This one controls “universal” settings – namely, titles,
   “super titles”, row and column labels, and a-b-c subplot labelling.

Now, instead of having to remember all of these verbose, one-liner
matplotlib commands like ``ax.set_title`` and ``ax.xaxis.tick_params``,
or even having to directly use verbose classes like the matplotlib
`~matplotlib.ticker` classes, `~proplot.axes.BaseAxes.format` lets
you change everything all at once. This basically eliminates the need
for boilerplate plotting code!

Also note the axes returned by `~proplot.subplots.subplots` function
are in a special `~proplot.subplots.axes_list` list. This lets you
call any method, including `~proplot.axes.BaseAxes.format`, on every
axes **simultaneously** (as in the below example).

.. code:: ipython3

    import proplot as plot
    f, axs = plot.subplots(ncols=2, nrows=2, share=False, span=False, tight=True, axwidth=1.5)
    axs.format(xlabel='x-axis', ylabel='y-axis', xlim=(0,10), xlocator=2,
              ylim=(0,4), ylocator=plot.arange(0,4), yticklabels=('a', 'bb', 'ccc', 'dd', 'e'),
              title='Axes title', titlepos='co', suptitle='Super title',
              abc=True, abcpos='il', abcformat='a.',
              ytickloc='both', yticklabelloc='both', ygridminor=True, xtickminor=False,
              collabels=['Column label 1', 'Column label 2'], rowlabels=['Row label 1', 'Row label 2'])



.. image:: showcase/showcase_20_0.svg


Automatic formatting
--------------------

ProPlot also adds back some of the convenience you get with the
`pandas` and `xarray` plotting commands. By default, when you plot a
`~pandas.DataFrame` or `~xarray.DataArray`, the x-axis label, y-axis
label, and title are updated from the metadata. To disable this
behavior, pass ``autoformat=False`` to `~proplot.subplots.subplots`.

The below examples showcase this feature for 1-dimensional and
2-dimensional datasets. It also demonstrates automatic generation of
colorbars and legends. For more on the ``colorbar`` and ``legend``
keywords, see `~proplot.wrappers.cmap_wrapper`,
`~proplot.wrappers.cycle_wrapper`, and the :ref:`Wrapped methods`
section. For more on panels, see the :ref:`Axes panels` section.

.. code:: ipython3

    import xarray as xr
    import numpy as np
    import pandas as pd
    import proplot as plot
    plot.rc['axes.formatter.timerotation']
    # DataArray
    # Must be column major since plot draws lines from columns of arrays
    data = np.sin(np.linspace(0, 2*np.pi, 20))[:,None] + np.random.rand(20,5).cumsum(axis=1)
    da = xr.DataArray(data, dims=('x', 'cat'), coords={
        'x':xr.DataArray(np.linspace(0,1,20), dims=('x',), attrs={'long_name':'distance', 'units':'km'}),
        'cat':xr.DataArray(np.arange(0,50,10), dims=('cat',), attrs={'long_name':'parameter', 'units':'K'})
        }, name='position series')
    # DataFrame
    ts = pd.date_range('1/1/2000', periods=20)
    data = (np.cos(np.linspace(0, 2*np.pi, 20))**4)[:,None] + np.random.rand(20,5)**2
    df = pd.DataFrame(data, index=ts, columns=['foo','bar','baz','foobar','barbaz'])
    df.name = 'time series'
    df.index.name = 'time (s)'
    df.columns.name = 'name'
    # Series
    series = pd.Series(np.random.rand(20).cumsum())
    # Figure
    f, axs = plot.subplots(ncols=2, share=False, span=False)
    axs.format(suptitle='Automatic subplot formatting')
    # Plot DataArray
    ax = axs[0]
    ax.plot(da, cycle=plot.shade('grass green', 0.3), lw=2, colorbar=True, colorbar_kw={'length':'2cm'})
    # Plot Dataframe
    ax = axs[1]
    ax.plot(df, cycle=plot.shade('caribbean green', 0.3), legend=True, legend_kw={'frameon':True}, lw=2)
    ax.format(xrotation=45)



.. image:: showcase/showcase_23_0.svg


.. code:: ipython3

    import xarray as xr
    import numpy as np
    import pandas as pd
    import proplot as plot
    from string import ascii_lowercase
    # DataArray
    data = 50*(np.sin(np.linspace(0, 2*np.pi, 20) + 0)**2) * np.cos(np.linspace(0, np.pi, 20)+np.pi/2)[:,None]**2
    da = xr.DataArray(data, dims=('plev','lat'), coords={
        'plev':xr.DataArray(np.linspace(1000,0,20), dims=('plev',), attrs={'long_name':'pressure', 'units':'hPa'}),
        'lat':xr.DataArray(np.linspace(-90,90,20), dims=('lat',), attrs={'units':'deg_N'}), # if long_name absent, variable name is used
        }, name='u', attrs={'long_name':'zonal wind', 'units':'m/s'})
    # DataFrame
    data = np.random.rand(20,20)
    df = pd.DataFrame(data.cumsum(axis=0).cumsum(axis=1), index=[*ascii_lowercase[:20]])
    df.name = 'funky data'
    df.index.name = 'index'
    df.columns.name = 'time (days)'
    # Figure
    # We must make room for the axes panels during subplots call!
    f, axs = plot.subplots(nrows=2, axcolorbars={1:'r', 2:'l'}, axwidth=2, share=False, span=False)
    axs.format(collabels=['Automatic subplot formatting']) # suptitle will look off center with the empty left panel
    # Plot DataArray
    ax = axs[0]
    ax.contourf(da, cmap='Tempo', colorbar='r')
    # Plot DataFrame
    ax = axs[1]
    ax.contourf(df, cmap='Speed', colorbar='l')
    ax.format(xtickminor=False)



.. image:: showcase/showcase_24_0.svg


Rc settings control
-------------------

A special object named `~proplot.rcmod.rc`, belonging to the
`~proplot.rcmod.rc_configurator` class, is created whenever you import
ProPlot. This object gives you advanced control over the look of your
plots. **Use** `~proplot.rcmod.rc` **as your one-stop shop for
changing global settings**.

To modify a setting for just one subplot, pass it to the
`~proplot.axes.BaseAxes.format` command. To reset everything to the
default state, use `~proplot.rcmod.rc_configurator.reset` (called by
default when a figure is rendered or saved by matplotlib; see
`~proplot.subplots.Figure`).

For more information, see the `~proplot.rcmod` documentation.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    # A bunch od different ways to update settings
    plot.rc.cycle = 'colorblind'
    plot.rc.linewidth = 1.5
    plot.rc.update({'fontname': 'DejaVu Sans'})
    plot.rc['figure.facecolor'] = 'w'
    plot.rc['axes.facecolor'] = 'gray5'
    # Make plot
    f, axs = plot.subplots(nrows=1, ncols=2, aspect=1, width=6,
                           span=0, wspace=0.5, sharey=2, hspace=0.7)
    N, M = 100, 6
    values = np.arange(1,M+1)
    for i,ax in enumerate(axs):
        data = np.cumsum(np.random.rand(N,M)-0.5, axis=0)
        lines = ax.plot(data, linewidth=3, cycle=('C0','C1',6)) # see "Changing the color cycle" for details
    axs.format(ytickloc='both', ycolor='blue7',
               hatch='xxx', hatchcolor='w',
               xlabel='x label', ylabel='y label',
               yticklabelloc='both',
               suptitle='Using "format" and "plot.rc" to apply new rc settings')
    ay = axs[-1].twinx()
    ay.format(ycolor='r', ylabel='secondary axis')
    ay.plot((np.random.rand(100)-0.2).cumsum(), color='r', lw=3)







.. image:: showcase/showcase_26_1.png
   :width: 540px
   :height: 260px


The `~proplot.rcmod.rc` object can also be used to change the default
font (see below). By default, ProPlot adds Helvetica and makes it the
new default. Helvetica is the MATLAB default – matplotlib normally does
not come packaged with Helvetica, but in my biased opinion it looks more
professional than the default “DejaVu Sans”. See the
`~proplot.fonttools` documentation for more info on fonts.

.. code:: ipython3

    import proplot as plot
    plot.rc.small =  8
    options = ['ultralight', 'light', 'normal', 'medium', 'demi', 'bold', 'extra bold', 'black']
    fonts = ['Helvetica', 'Helvetica Neue', 'DejaVu Sans', 'Bitstream Vera Sans', 'Verdana', 'Tahoma', 'Arial', 'Geneva']
    f, axs = plot.subplots(ncols=4, nrows=len(fonts)//4, share=False, span=False,
                           axwidth=1.5, axheight=2, wspace=0.5, hspace=0.5)
    #, 'Times New Roman', 'Palatino', 'Inconsolata', 'Myriad Pro']
    for ax,font in zip(axs,fonts):
        plot.rc['fontname'] = font
        math  = r'$\alpha\beta + \gamma\delta \times \epsilon\zeta \cdot \eta\theta$'
        math += ('\n' + r'$\Sigma\kappa\lambda\mu\pi\rho\sigma\tau\psi\phi\omega$')
        ax.text(0.5, 0, math + '\n' + 'The quick brown fox\njumps over the lazy dog.\n0123456789\n!@#$%^&*()[]{};:,./?',
                weight='normal', ha='center', va='bottom')
        ax.format(xlabel='xlabel', ylabel='ylabel', suptitle='Table of font names')
        for i,option in enumerate(options):
            if option in ('italic', 'oblique'):
                kw = {'style':option, 'weight':'normal'} # otherwise defaults to *lightest* one!
            elif option in ('small-caps',):
                kw = {'variant':option}
            else:
                kw = {'weight':option}
            kw.update({'stretch':'normal'})
            ax.text(0.03, 0.97 - (i*1.2*(plot.rc['small']/72)/ax.height), f'{option}', ha='left', va='top', **kw)
            ax.text(0.97, 0.97 - (i*1.2*(plot.rc['small']/72)/ax.height), f'{font[:14].strip()}',   ha='right', va='top', **kw)



.. image:: showcase/showcase_28_0.svg


