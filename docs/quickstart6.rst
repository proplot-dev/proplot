Plotting wrappers
=================

New features have been added to various matplotlib plotting commands
thanks to a set of wrapper functions. These features are a strict
*superset* of the existing matplotlib API – if you want, you can use
plotting commands exactly as you always have. This section documents
these wrapper functions. For details, see the `~proplot.axes`
documentation.

You should also see :ref:`Making your own colormaps` and
:ref:`Making your own color cycles`, which explain how
`~proplot.wrappers.cmap_wrapper` and
`~proplot.wrappers.cycle_wrapper` can be used to create and apply new
colormaps and property cyclers on-the-fly.

Colormap normalizers
--------------------

`~proplot.wrappers.cmap_wrapper` assigns the
`~proplot.colortools.BinNorm` “meta-normalizer” as the data normalizer
for all plotting commands involving colormaps. This permits discrete
``levels`` even for commands like `~matplotlib.axes.Axes.pcolor` and
`~matplotlib.axes.Axes.pcolormesh`. `~proplot.colortools.BinNorm`
also ensures that colorbar colors span the entire colormap range,
independent of the ``extend`` setting, and that color levels on the ends
of colorbars for “cyclic” colormaps are distinct.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    f, axs = plot.subplots(ncols=2, axwidth=1.5, axcolorbars={1:'l', 2:'r'})
    cmap = 'orange5'
    data = np.random.rand(15,15)
    axs.format(suptitle='Pcolor with levels demo')
    ax = axs[0]
    ax.pcolor(data, cmap=cmap, colorbar='l', vmin=0, vmax=1, levels=200, colorbar_kw={'ticks':0.2})
    ax.format(title='Ambiguous values', yformatter='null')
    ax = axs[1]
    ax.pcolor(data, cmap=cmap, colorbar='r', levels=np.linspace(0,1,6), colorbar_kw={'ticks':0.2})
    ax.format(title='Discernible values')



.. image:: quickstart/quickstart_143_0.svg


.. code:: ipython3

    import proplot as plot
    import numpy as np
    f, axs = plot.subplots(ncols=5, width=8, wratios=(5,3,3,3,3), axcolorbars='b')
    axs.format(suptitle='Demo of colorbar color-range standardization')
    levels = plot.arange(0,360,45)
    data = (20*(np.random.rand(20,20) - 0.4).cumsum(axis=0).cumsum(axis=1)) % 360
    ax = axs[0]
    ax.pcolormesh(data, levels=levels, cmap='phase', extend='neither', colorbar='b')
    ax.format(title='Cyclic map with separate ends')
    for ax,extend in zip(axs[1:], ('min','max','neither','both')):
        ax.pcolormesh(data, levels=levels, cmap='spectral', extend=extend, colorbar='b', colorbar_kw={'locator':90})
        ax.format(title=f'Map with extend={extend}')



.. image:: quickstart/quickstart_144_0.svg


If you pass unevenly spaced ``levels``, the
`~proplot.colortools.LinearSegmentedNorm` normalizer is applied by
default. This results in even color gradations across *indices* of the
level list, no matter their spacing. To use an arbitrary colormap
normalizer, just pass ``norm`` and optionally ``norm_kw`` to a command
wrapped by `~proplot.wrappers.cmap_wrapper`. These arguments are
passed to the `~proplot.colortools.Norm` constructor.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    f, axs = plot.subplots(axcolorbars='b', ncols=2, axwidth=2.5, aspect=1.5)
    data = 10**(2*np.random.rand(20,20).cumsum(axis=0)/7)
    ticks = [5, 10, 20, 50, 100, 200, 500, 1000]
    for i,(norm,title) in enumerate(zip(('linear','segments'),('Linear normalizer','LinearSegmentedNorm (default)'))):
        m = axs[i].contourf(data, levels=ticks, extend='both', cmap='Mako', norm=norm, colorbar='b')
        axs[i].format(title=title)
    axs.format(suptitle='Level normalizers demo')



.. image:: quickstart/quickstart_146_0.svg


Finally, there is a new `~proplot.colortools.MidpointNorm` class that
warps your colormap so that its midpoint lies on some central data
value, no matter the minimum and maximum colormap colors. Again, to use
an arbitrary colormap normalizer, just pass ``norm`` and optionally
``norm_kw`` to a command wrapped by `~proplot.wrappers.cmap_wrapper`.
These arguments are passed to the `~proplot.colortools.Norm`
constructor.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    data1 = (np.random.rand(20,20) - 0.43).cumsum(axis=0)
    data2 = (np.random.rand(20,20) - 0.57).cumsum(axis=0)
    f, axs = plot.subplots(ncols=2, axwidth=2.5, aspect=1.5, axcolorbars='b')
    cmap = plot.Colormap('Moisture', cut=0.1)
    axs.format(suptitle='Midpoint normalizer demo')
    axs[0].contourf(data1, norm='midpoint', cmap=cmap, colorbar='b')
    axs[0].format(title='Skewed positive data')
    axs[1].contourf(data2, norm='midpoint', cmap=cmap, colorbar='b')
    axs[1].format(title='Skewed negative data')



.. image:: quickstart/quickstart_148_0.svg


Heatmaps and labeling
---------------------

The new `~proplot.axes.BaseAxes.heatmap` command calls
`~matplotlib.axes.Axes.pcolormesh` and applies default formatting that
is suitable for heatmaps: no minor ticks, no gridlines, and major ticks
at the center of each box.

You can also add labels to `~matplotlib.axes.Axes.pcolor`,
`~matplotlib.axes.Axes.pcolormesh`,
`~proplot.axes.BaseAxes.heatmap`, `~matplotlib.axes.Axes.contour`,
and `~matplotlib.axes.Axes.contourf` plots, thanks to
`~proplot.wrappers.cmap_wrapper`. Just pass the ``labels=True``
keyword argument, and ProPlot will draw contour labels with
`~matplotlib.axes.Axes.clabel` or grid box labels with
`~matplotlib.axes.Axes.text`. The label format can be changed by
passing a ``labels_kw`` dictionary of settings (e.g.
``labels_kw={'fontsize':12}``) and with the ``precision`` keyword arg.
Label colors are automatically chosen based on the luminance of the
underlying box or contour color.

.. code:: ipython3

    import proplot as plot
    import pandas as pd
    import numpy as np
    f, axs = plot.subplots(axwidth=2, ncols=2, span=False, share=False)
    data = np.random.rand(6,6)
    data = pd.DataFrame(data, index=pd.Index(['a','b','c','d','e','f']))
    axs.format(suptitle='Labels demo')
    ax = axs[0]
    m = ax.heatmap(data, cmap='rocket', labels=True, precision=2, labels_kw={'weight':'bold'})
    ax.format(xlabel='xlabel', ylabel='ylabel', title='Heatmap plot with bold labels')
    ax = axs[1]
    m = ax.contourf(data.cumsum(axis=0), labels=True, cmap='rocket', labels_kw={'weight':'bold'})
    ax.format(xlabel='xlabel', ylabel='ylabel', title='Contourf plot with bold labels')



.. image:: quickstart/quickstart_151_0.svg


Easy error bars
---------------

Thanks to the `~proplot.wrappers.add_errorbars` wrapper, you can now
add error bars when using the `~matplotlib.axes.Axes.plot`,
`~matplotlib.axes.Axes.scatter`, `~matplotlib.axes.Axes.bar`,
`~matplotlib.axes.Axes.barh`, and `~matplotlib.axes.Axes.violinplot`
methods. If you pass 2D arrays of data to these commands with
``means=True`` or ``medians=True``, the *means or medians* of each
column are drawn as points, lines, or bars, and error bars represent the
*spread* in each column. You can draw both thin “bars” with optional
whiskers, and thick “boxes” overlayed on top of these bars. You can also
pass error bar coordinates manually with the ``bardata`` and ``boxdata``
keyword args. See `~proplot.wrappers.add_errorbars` for details.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    import pandas as pd
    plot.rc['title.loc'] = 'uc'
    plot.rc['axes.ymargin'] = plot.rc['axes.xmargin'] = 0.05
    data = np.random.rand(20,8).cumsum(axis=0).cumsum(axis=1)[:,::-1] + 20*np.random.normal(size=(20,8)) + 30
    f, axs = plot.subplots(nrows=3, aspect=1.5, axwidth=3, span=False, share=False, hratios=(2,1,1))
    axs.format(suptitle='Error bars with various plotting commands')
    # Asking add_errorbars to calculate bars
    ax = axs[0]
    obj = ax.barh(data, color='red orange', means=True)
    ax.format(title='Column statistics')
    # Showing a standard deviation range instead of percentile range
    ax = axs[1]
    ax.scatter(data, color='goldenrod', marker='x', markersize=50, barcolor='gray8',
               medians=True, barstd=True, barrange=(-1,1), barzorder=0, boxes=False, capsize=2)
    # Supplying error bar data manually
    ax = axs[2]
    boxdata = np.percentile(data, (25,75), axis=0)
    bardata = np.percentile(data, (5,95), axis=0)
    ax.plot(data.mean(axis=0), lw=2, barlw=1, boxmarker=False, edgecolor='gray6', color='yellow orange',
            boxdata=boxdata, bardata=bardata)
    # Formatting
    axs[0].format(ylabel='column number', title='Bar plot')
    axs[1].format(title='Scatter plot')ss
    axs[2].format(title='Line plot')
    axs[1:].format(xlabel='column number', xticks=1)



.. image:: quickstart/quickstart_154_0.svg


Parametric plots
----------------

`~matplotlib.axes.Axes.plot` now accepts a ``cmap`` keyword – this
lets you draw line collections that map individual segments of the line
to individual colors. This can be useful for drawing “parametric” plots,
where you want to indicate the time or some other coordinate at each
point on the line. See `~proplot.axes.BaseAxes.cmapline` for details.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    f, axs = plot.subplots(span=False, share=False, ncols=2, wratios=(2,1), axcolorbars='b', axwidth='5cm', aspect=(2,1))
    ax = axs[0]
    m = ax.plot((np.random.rand(50)-0.5).cumsum(), np.random.rand(50),
                cmap='thermal', values=np.arange(50), lw=7, extend='both')
    ax.format(xlabel='xlabel', ylabel='ylabel', title='Line with smooth color gradations', titleweight='bold')
    ax.colorbar(m, loc='b', label='parametric coordinate', locator=5)
    N = 12
    ax = axs[1]
    values = np.arange(1, N+1)
    radii = np.linspace(1,0.2,N)
    angles = np.linspace(0,4*np.pi,N)
    x = radii*np.cos(1.4*angles)
    y = radii*np.sin(1.4*angles)
    m = ax.plot(x, y, values=values,
                linewidth=15, interp=False, cmap='thermal')
    ax.format(xlim=(-1,1), ylim=(-1,1), title='With step gradations', titleweight='bold',
              xlabel='cosine angle', ylabel='sine angle')
    ax.colorbar(m, loc='b', locator=None, label=f'parametric coordinate')







.. image:: quickstart/quickstart_157_1.svg


Area plots
----------

Make area plots with the convenient aliases
`~proplot.axes.BaseAxes.area` and `~proplot.axes.BaseAxes.areax`.
These point to the `~matplotlib.axes.Axes.fill_between` and
`~matplotlib.axes.Axes.fill_betweenx` methods, which are wrapped with
`~proplot.wrappers.fill_between_wrapper` and
`~proplot.wrappers.fill_betweenx_wrapper`.

The wrappers enable “stacking” successive columns of a 2D input array
like in `pandas`. They also add a new “``negpos``” keyword for
creating area plots that change color when the fill boundaries cross
each other. The most common use case for this is highlighting negative
and positive area underneath a line, as shown below.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    f, axs = plot.subplots(array=[[1,2],[3,3]], hratios=(1,0.8), span=False, share=0)
    axs.format(xlabel='xlabel', ylabel='ylabel', suptitle='Area plot demo')
    data = np.random.rand(5,3).cumsum(axis=0)
    ax = axs[0]
    ax.areax(np.arange(5), data, data + np.random.rand(5)[:,None], alpha=0.5,
            legend='uc', legend_kw={'center':True, 'ncols':2, 'labels':['z','y','qqqq']},
            )
    ax.format(title='Fill between columns')
    ax = axs[1]
    ax.area(np.arange(5), data, stacked=True, alpha=0.8,
            legend='ul', legend_kw={'center':True, 'ncols':2, 'labels':['z','y','qqqq']},
            )
    ax.format(title='Stack between columns')
    ax = axs[2]
    data = 5*(np.random.rand(20)-0.5)
    ax.area(data, negpos=True, negcolor='blue7', poscolor='red7')
    ax.format(title='Negative and positive data', xlabel='xlabel', ylabel='ylabel')



.. image:: quickstart/quickstart_160_0.svg


Bar plots
---------

`~proplot.wrappers.bar_wrapper` and
`~proplot.wrappers.cycle_wrapper` make it easier to generate useful
bar plots. You can now pass 2D arrays to `~matplotlib.axes.Axes.bar`
or `~matplotlib.axes.Axes.barh`, and columns of data will be grouped
or stacked together. And if *x* coordinates are not provided, default
coordinates are applied, just like with `~matplotlib.axes.Axes.plot`.
See `~proplot.wrappers.bar_wrapper` for details.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    import pandas as pd
    plot.rc.titleloc = 'uc'
    plot.rc.margin = 0.05
    f, axs = plot.subplots(nrows=2, aspect=2, axwidth=3, span=False, share=False)
    data = np.random.rand(5,5).cumsum(axis=0).cumsum(axis=1)[:,::-1]
    data = pd.DataFrame(data, columns=pd.Index(np.arange(1,6), name='column'), index=pd.Index(['a','b','c','d','e'], name='row idx'))
    ax = axs[0]
    obj = ax.bar(data, cycle='Reds', cycle_kw={'left':0.2}, colorbar='ul', colorbar_kw={'frameon':False})
    ax.format(xlocator=1, xminorlocator=0.5, ytickminor=False, title='Side-by-side', suptitle='Bar plot wrapper demo')
    ax = axs[1]
    obj = ax.barh(data.iloc[::-1,:], cycle='Grays', legend='ur', stacked=True)
    ax.format(title='Stacked')



.. image:: quickstart/quickstart_163_0.svg


Box plots and violins
---------------------

`~matplotlib.axes.Axes.boxplot` and
`~matplotlib.axes.Axes.violinplot` are now wrapped with
`~proplot.wrappers.boxplot_wrapper`,
`~proplot.wrappers.violinplot_wrapper`, and
`~proplot.wrappers.cycle_wrapper`, making it much easier to plot
distributions of data with aesthetically pleasing default settings and
automatic axis labeling.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    import pandas as pd
    f, axs = plot.subplots(ncols=2)
    data = np.random.normal(size=(20,5)) + 2*(np.random.rand(20,5)-0.5)
    data = pd.DataFrame(data, columns=pd.Index(['a','b','c','d','e'], name='xlabel'))
    ax = axs[0]
    obj1 = ax.boxplot(data, lw=0.7, marker='x', fillcolor='gray5', medianlw=1, mediancolor='k')#, boxprops={'color':'C0'})#, labels=data.columns)
    ax.format(title='Box plots', titleloc='uc')
    ax = axs[1]
    obj2 = ax.violinplot(data, lw=0.7, fillcolor='gray7', means=True)
    ax.format(title='Violin plots', titleloc='uc')
    axs.format(ymargin=0.1, xmargin=0.1, suptitle='Boxes and violins demo')



.. image:: quickstart/quickstart_166_0.svg


Scatter plots
-------------

Thanks to `~proplot.wrappers.scatter_wrapper` and
`~proplot.wrappers.cycle_wrapper`, `~matplotlib.axes.Axes.scatter`
now accepts 2D arrays, just like `~matplotlib.axes.Axes.plot`, and
successive calls to `~matplotlib.axes.Axes.scatter` can apply property
cycle keys other than ``color`` – for example, ``marker`` and
``markersize``. `~matplotlib.axes.Axes.scatter` also now optionally
accepts keywords that look like the `~matplotlib.axes.Axes.plot`
keywords, which is a bit less confusing. You can also pass colormaps to
`~matplotlib.axes.Axes.scatter` just as with matplotlib.

.. code:: ipython3

    import proplot as plot
    import numpy as np
    import pandas as pd
    plot.rc.reset()
    f, axs = plot.subplots(ncols=2, share=1)
    x = (np.random.rand(20)-0).cumsum()
    data = (np.random.rand(20,4)-0.5).cumsum(axis=0)
    data = pd.DataFrame(data, columns=pd.Index(['a','b','c','d'], name='label'))
    # Scatter demo
    ax = axs[0]
    ax.format(title='New prop cycle properties', suptitle='Scatter plot demo')
    obj = ax.scatter(x, data, legend='ul', cycle='538', legend_kw={'ncols':2},
                    cycle_kw={'marker':['x','o','x','o'], 'markersize':[5,10,20,30]})
    ax = axs[1]
    ax.format(title='Scatter colormap with colorbar')
    data = (np.random.rand(2,100)-0.5)
    obj = ax.scatter(*data, color=data.sum(axis=0), size=10*(data.sum(axis=0)+1),
                     marker='*', cmap='fire', colorbar='ll', colorbar_kw={'locator':0.5, 'label':'label'})
    axs.format(xlabel='xlabel', ylabel='ylabel')



.. image:: quickstart/quickstart_169_0.svg
