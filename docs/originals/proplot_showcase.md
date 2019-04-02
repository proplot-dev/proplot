<!-- See: https://predictablynoisy.com/jekyll-markdown-nbconvert -->
<!-- This is a simple template for trimming everything except cells -->
<!-- and the figures they spit out. All remaining output is stripped. -->
<!-- 
This is provided as a *reference*, to see what the 'nbtemplate' template
points to, and see the available block names that can be modified. Apparently
the super() block says to use the parent 'extends' template method.
-->

# The basics
Use `plot.subplots()` to generate plots of arbitrary complexity, assign map
projections to particular axes, add panels to the edges of subplots or the
entire figure, and designate plot dimensions compatible for particular
publications. Note that figures drawn in *interactive sessions* will have gray
backgrounds with white axes, while *saved* figures will be transparent (disable
this by passing `transparent=False` to `fig.save()`).

Every axes generated with `plot.subplots()` is a special subclass of the
`matplotlib.axes.Axes` class, with several new methods introduced. The most
important of these is the `format()` method. This command is extremely powerful,
and can be used to create highly customized figures -- it is best demonstrated
by example (see below).

## Share and span
Pass the `nrows` and `ncols` arguments to `subplots()` to create simple subplot
grids.

This library expands the builtin "[shared
axis](https://matplotlib.org/examples/pylab_examples/shared_axis_demo.html)"
matplotlib feature, using `subplots()` arguments `sharex`, `sharey`, and `share`
(both axes). Use `share[x|y]=0` for no axis sharing `share[x|y]=1` for sharing
axis labels, but not limits or tick labels, `share[x|y]=2` for sharing axis
labels and limits, but not tick labels, and `share[x|y]=3` for sharing axis
labels, limits, and tick labels. There is also a new feature for making
"subplot-spanning" x and y labels, toggled with the `spanx`, `spany`, or `span`
arguments.

A worked example is below, where the **y-axes** are "shared" and the **x-axes**
are labelled with a "spanning" axis label.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
N = 50
M = 40
colors = plot.colors('grays_r', M, x=(0.1, 0.8))
for share in (0,1,2,3):
    f, axs = plot.subplots(ncols=4, aspect=1.2, wspace=0.5, axwidth=1.2, sharey=share, spanx=share//2)
    gen = lambda scale: scale*(np.random.rand(N,M)-0.5).cumsum(axis=0)[N//2:,:]
    for ax,scale,color in zip(axs,(1,3,7,0.2),('gray9','gray7','gray5','gray3')):
        array = gen(scale)
        for l in range(array.shape[1]):
            ax.plot(array[:,l], color=colors[l])
        ax.format(suptitle=f'Axis-sharing level: {share}, spanning labels {["off","on"][share//2]}', ylabel='y-label', xlabel='x-axis label')
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Adjusting gridspec.
Adjusting gridspec.
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_2_1.svg)



![svg](/tools/files/proplot_showcase_files/proplot_showcase_2_2.svg)



![svg](/tools/files/proplot_showcase_files/proplot_showcase_2_3.svg)



![svg](/tools/files/proplot_showcase_files/proplot_showcase_2_4.svg)



```python
import proplot as plot
import numpy as np
plot.nbsetup()
plot.rc.cycle = 'Set4'
titles = ['With redundant labels', 'Without redundant labels']
for mode in (0,1):
    f, axs = plot.subplots(nrows=4, ncols=4, share=3*mode, span=1*mode, axwidth=1,
                           wspace=0.2 + 0.4*(1-mode), hspace=0.15 + 0.25*(1-mode))
    for ax in axs:
        ax.plot((np.random.rand(100,20)-0.4).cumsum(axis=0))
    axs.format(xlabel='x-label', ylabel='y-label', suptitle=titles[mode], abc=mode)
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_3_1.svg)



![svg](/tools/files/proplot_showcase_files/proplot_showcase_3_2.svg)


## Subplot arrays
Set up an arbitrarily complex grid of subplots using a 2D array of integers (or
iterable of iterables). Use `0` for empty spaces, and `1` to `N` for unique
subplots, `N` being the number of subplots you want. You can think of this array
as a "picture" of the grid you want. The below example demonstrates this nicely.

The list of `Axes` instances returned by `subplots()` is actually a **special
class** called `axes_list`, ordered by the numbering you used in the 2D array
(if you used the `nrows` or `ncols` arguments instead, default numbering is row-
major order). You can bulk-call any method across several axes by accessing that
method attribute on the `axes_list` -- this is done with the `format()` method
below. You can also use `axs.item` to retrieve a *list* of `item` attributes
from each axes in the `axes_list`.


```python
# Arbitrarily complex array of subplots, with shared/spanning x/y axes detected automatically
import proplot as plot
import numpy as np
plot.nbsetup()
f, axs = plot.subplots([[1, 1, 2], [1, 1, 6], [3, 4, 4], [3, 5, 5]],
                       span=1, share=3,
                       wspace=0.6, hspace=0.5, width=5)
axs.format(suptitle='Complex subplot grid', xlabel='time (seconds)', ylabel='temperature (K)', abc=True)
axs[0].plot(2*(np.random.rand(100,5)-0.5).cumsum(axis=0))
```

<!--
Configured ipython notebook.

```
{:.output_stream}
```
-->




<!--
{:.output_data_text}
```
[<matplotlib.lines.Line2D at 0xb1a5ccc88>,
 <matplotlib.lines.Line2D at 0xb19d38198>,
 <matplotlib.lines.Line2D at 0xb19d38240>,
 <matplotlib.lines.Line2D at 0xb19d386d8>,
 <matplotlib.lines.Line2D at 0xb19d38668>]
```
-->



<!--
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_5_3.svg)


## Outer panels, formatting
Easily label rows/columns of your plot, add figure titles, add colorbars/legends
to the perimeter of the subplot region, label subplots with a-b-c enumeration,
and reposition titles. In this example, we use outer "panels" to draw colorbars.
Each panel may *span arbitrary contiguous rows and columns*. For more
information on panels, see the
[documentation](https://lukelbd.github.io/tools/proplot/doc).


```python
# Multiple subplots, long axes
import proplot as plot
import numpy as np
plot.nbsetup()
f, axs = plot.subplots(tight=True, spany=False, sharey=3, sharex=1,
                       nrows=3, ncols=3, axwidth=1.5, aspect=1,
                       wratios=[2,1,1], hratios=[2,1,1],
                       bottom=0.5, left=0.5, lspace=0.5,
                       hspace=0.3, wspace=(0.2, 0.6),
                       # hspace=(0.1, 0.4), wspace=(0.1, 0.4),
                       bottompanel=True, rightpanels=[1,2,2])
m = axs[0].contourf(np.random.rand(10,10).cumsum(axis=0), rowmajor=True, extend='both')
# axs[:3].format(title='Minor titles')
axs.format(abc=True, abcpos='li', abcformat='a.',
           suptitle='SuperTitle is automatically offset and centered above main axes',
           title='Inner title', titlepos='inside', # title_kw={'fancy':True},
           collabels=['Column A', 'Column B', 'Column C'], collabels_kw=dict(color='k', weight='bold'),
           rowlabels=['Row 1', 'Row 2', 'Row 3'], rowlabels_kw=dict(color='k', weight='bold'),
           xlabel='xlabel', ylabel='ylabel')
# axs[-1].format(color='r', linewidth=1.1)
axs[-1].format(linewidth=1.1, color='r')
f.bottompanel.colorbar(m, length=0.9, cgrid=True, cformatter='none', clocator='none')
res = f.rightpanel[:2].colorbar(m, clabel='clabel', ctickminor=True, clocator=1, cminorlocator=0.5, extend='neither') # draws two colorbars simultaneously
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_7_1.svg)


## Inner panels, rc settings
Modify global settings with `plot.rc['prop']` or `plot.rc.prop`. This includes
`rcParams` settings (i.e. builtin matplotlib global settings), custom
`rcSpecial` settings, and some bulk `rcGlobals` settings that apply to multiple
other settings. See the
[documentation](https://lukelbd.github.io/tools/proplot/doc) for more
information settings configuration in ProPlot.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
plot.rc.linewidth = 1.2
f, axs = plot.subplots(nrows=1, ncols=2, aspect=0.8, width=6,
                       spanx=1, spany=0, right=0.6, wspace=0.5,
                       sharex=0, sharey=2, hspace=0.7, bottom=0.5,
                       innerpanels='b', bottomcolorbar=True,
                      )
                     # innerpanels=True, whichpanels='b')
N, M = 100, 6
values = np.arange(1,M+1)
for i,ax in enumerate(axs):
    plot.rc.cycle = ['C0','C1',6]
    data = np.cumsum(np.random.rand(N,M)-0.5, axis=0)
    lines = ax.plot(data, linewidth=2)
    ax.bottompanel.plot(data.mean(axis=1), color='gray7', lw=2)
axs.format(ytickloc='both', ycolor='blue7', xlabel='spanning x label', ylabel='ylabel', abc=True, abcpos='il',
           yticklabelloc='both',
           suptitle='Various features demonstrated below')
ay = axs[-1].twinx()
ay.format(ycolor='r', ylabel='secondary axis')
ay.plot((np.random.rand(100)-0.2).cumsum(), color='r', lw=2)
f.bottompanel.colorbar(lines, values=values, length=0.7, extend='both', clocator=values, clabel='time series no.')
```

<!--
Configured ipython notebook.

```
{:.output_stream}
```
-->




<!--
{:.output_data_text}
```
(<matplotlib.axes._subplots.AxesSubplot at 0xb1d620518>,
 <matplotlib.colorbar.Colorbar at 0xb1d6175f8>)
```
-->



<!--
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_9_3.svg)


## Font selection
Easily switch between different fonts using the `fontname` rc property. The
`ttf` files from several fonts are distributed with this package, and can be
locally installed to your matplotlib distribution using `plot.install_fonts()`
(may require restarting iPython session). Note `plot.fonts` is a list of
available fonts, `plot.fonts_os` is a list of OS-provided fonts, and
`plot.fonts_mpl` is a list of fonts packaged with matplotlib (or added with
`install_fonts`).


```python
import proplot as plot
plot.nbsetup()
plot.rc['small'] = plot.rc['large'] = 10
plot.rc['fontname'] = 'Helvetica'
f, axs = plot.subplots(ncols=4, nrows=3, share=False, span=False,
                       axwidth=2.0, aspect=0.85, wspace=0.5, hspace=0.5)
# options = ['ultralight', 'light', 'normal', 'regular', 'book', 'medium', 'roman',
#            'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black',
#            'italic', 'oblique'] # remove redundancies below
options = ['ultralight', 'light', 'normal', 'medium', 'demi', 'bold', 'extra bold', 'black']
fonts = ['Helvetica', 'Helvetica Neue', 'DejaVu Sans', 'Bitstream Vera Sans', 'Verdana', 'Tahoma',
         'Arial', 'Geneva', 'Times New Roman', 'Palatino', 'Inconsolata', 'Myriad Pro'] #Comic Sans MS', 'Myriad Pro']
for ax,font in zip(axs,fonts):
    plot.rc['fontname'] = font
    math  = r'$\alpha\beta + \gamma\delta \times \epsilon\zeta \cdot \eta\theta$'
    math += ('\n' + r'$\Sigma\kappa\lambda\mu\pi\rho\sigma\tau\psi\phi\omega$')
    ax.text(0.5, 0, math + '\n' + 'The quick brown fox\njumps over the lazy dog.\n0123456789\n!@#$%^&*()[]{};:,./?',
            weight='normal', ha='center', va='bottom')
    ax.format(xlabel='xlabel', ylabel='ylabel')#, title=font, titlepos='il', title_kw={'border':False, 'weight':'bold'}) #, rc_kw={'fontname':font})
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
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_11_1.svg)


## Cartopy vs. Basemap
Here we can see how Cartopy's low-level integration with the matplotlib API
shines. With basemap, your data must simply be transformed to map-projection
coordinates. With cartopy, the underlying plotting tools operate in map-
projection coordinates.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
# First make figure
f, axs = plot.subplots(ncols=2, nrows=2, width=7, hspace=0.2, wspace=0.3, top=0.5,
                       bottomcolorbars=True, bwidth=0.2, bottom=0.2,
                       proj='hammer', proj_kw={'lon_0':0},
                       # basemap=False,
                       basemap={(1,3):False, (2,4):True},
                       )
offset = 20
x = plot.arange(-180+offset,180+offset-1,60)
y = plot.arange(-60,60+1,30)
data = np.random.rand(len(x), len(y))
for ax,p,pcolor,basemap in zip(axs,range(4),[1,1,0,0],[0,1,0,1]):
    # adfdas
    m = None
    cmap = ['sunset', 'sunrise'][basemap]
    levels = [0, .3, .5, .7, .9, 1]
    levels = np.linspace(0,1,11)
    if pcolor:
        m = ax.pcolorpoly(x, y, data, levels=levels, cmap=cmap, extend='both', extremes=True)
        ax.scatter(np.random.rand(5,5)*180, 180*np.random.rand(5,5))
    if not pcolor:
        m = ax.contourf(x, y, data, levels=levels, cmap=cmap, extend='both', extremes=False)
        ax.scatter(np.random.rand(5,5)*180, 180*np.random.rand(5,5))
    ax.format(facecolor='gray2', suptitle='Hammer projection in different mapping frameworks', collabels=['Cartopy', 'Basemap'])
    if p<2:
        ax, c = f.bottompanel[p].colorbar(m, clabel='values', ctickminor=False)
    # print(p, ax._sharex, ax._sharey, list(ax._shared_x_axes))
    # if p==2:
        # raise Exception
```

<!--
Configured ipython notebook.
Warning: Cannot label meridians on Hammer basemapWarning: Cannot label meridians on Hammer basemapResetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_13_1.svg)


Meridian and parallel labelling only work for stereographic and Mercator
projections in cartopy. However, the cartopy API is much more flexible and much
more powerful. Even complex plotting algorithms like `tricontour` work with
cartopy. Another custom colormap is used below.


```python
import proplot as plot
plot.nbsetup()
import numpy as np
f, axs = plot.subplots(ncols=2, width=7, proj={1:'merc', 2:'nplaea'},
                       wspace=0.5, basemap={1:False, 2:True},
                       proj_kw={1:{'lon_0':0}, 2:{'lon_0':0, 'boundinglat':5}}, left=0.4, right=0.4, bottom=0.2)
# First the tricolor cartopy plot
axs.set_adjustable('box')
ax = axs[0]
np.random.seed(3498)
x, y = np.random.uniform(size=(100, 2)).T
z = np.exp(-x**2 - y**2)
x = (x-0.5)*360
y = (y-0.5)*180
levels = np.linspace(0, 1, 100)
cnt = ax.tripcolor(x, y, z, levels=levels, cmap='Sea')
ax.format(title='Tricontour plot', xlabels='b', xlocator=60, ylocator=20)
# Next the basemap one
ax = axs[1]
N = 20
x = np.linspace(-180, 180, N)
x = x[:-1] # smooth transition across cutoff
y = np.linspace(-70, 70, N)
levels = np.linspace(0, 1, 100)
ax.format(title='Basemap plot', xlocator=plot.arange(-180,180,60), ylocator=plot.arange(-80,80,20),
          lonlabels='lrb', latlabels='')
cnt = ax.contourf(x, y, np.random.rand(len(x), len(y)).cumsum(axis=0), cmap='Sea', levels=20)
```

<!--
Configured ipython notebook.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_15_1.svg)


# New axis scales
This package also provides some special axis "scales", along with a tool for
creating arbitrary scales with "jumps" and "zooms".

## Latitude scales
The `sine` scale creates a geographically "area-weighted" latitude axis. The
`mercator` scale creates an axis in Mercator latitude coordinates, which is
occasionally useful [scientific
contexts](https://journals.ametsoc.org/doi/full/10.1175/JAS-D-11-039.1).


```python
import proplot as plot
import numpy as np
plot.nbsetup()
plot.rc.update(color='gray7', facehatch='xxxx')
f, axs = plot.subplots(ncols=2, width=7, share=0, span=0, wspace=0.7, left=0.6)
n = 30
x = np.linspace(-180,180,n)
y = np.linspace(-85,85,n) # note sine just truncated values not in [-90,90], but Mercator transformation can reflect them
y2 = np.linspace(-85,85,n) # for pcolor
for i,(ax,scale,color) in enumerate(zip(axs,['mercator','sine'],['sky blue','coral'])):
    ax = axs[i-1]
    ax.plot(x, y, '-', color=color, lw=4)
    data = np.random.rand(len(x), len(y2))
    ax.pcolormesh(x, y2, data, cmap='grays', cmap_kw={'right': 0.8}) # use 'right' to trim the colormap from 0-1 color range to 0-0.8 color range
    ax.format(xlabel='longitude', ylabel='latitude', title=scale.title() + '-latitude y-axis', yscale=scale,
              ytickloc='left', suptitle='Projection coordinate y-axes',
              xformatter='deglon', yformatter='deglat', grid=False,
              xscale='linear', xlim=None, ylim=(-85,85))
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_18_1.svg)


## Inverse scale
A scale useful primarily where you'd like to show the wavenumber and wavelength
on the same axis.


```python
# Plot the response function for an imaginary 5-day lowpass filter
import proplot as plot
import numpy as np
plot.nbsetup()
plot.rc['axes.ymargin'] = 0
cutoff = 0.3
x = np.linspace(0.01,0.5,1000) # in wavenumber days
response = (np.tanh(-((x - cutoff)/0.03)) + 1)/2 # imgarinary response function
f, ax = plot.subplots(aspect=(3,1), width=6)#, tight=False, top=2)
ax.fill_between(x, 0, response, hatch='xxx', facecolor='none', edgecolor='gray8', lw=1, clip_on=True)
ax.axvline(cutoff, lw=2, dashes=(0.2,2), color='red')
ax.format(xlabel='wavenumber (days$^{-1}$)', ylabel='response', grid=False)
axy = ax.twiny()
axy.format(xlim=(1/max(x), 1/min(x)), xlocator=np.array([20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05]),
          xscale='inverse', xlabel='period (days)',
          title='Title automatically offset above axis labels', titlepos='oc',
          suptitle='SuperTitle above everything', 
          )
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_20_1.svg)


## Pressure and height scales
Scale a height coordinate to correspond linearly with pressure using
`[x|y]scale='height'`. Scale a pressure coordinate to correspond linearly with
height using `[x|y]scale='pressure'`. Note the scale height assumed for these
conversions is 7km -- change this by using `[x|y]scale=('height', scale_height)`
or `[x|y]scale=('pressure', scale_height)`.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
cutoff = 0.1
f, axs = plot.subplots(aspect=(1,2.5), ncols=2,
                       bottom=0.4,
                       span=False, share=False, wspace=1, width=5, bottomlegend=True)
N = 500
H = 7.0
p0 = 1000.0
ylim = np.array([0, 25])
ylims = [ylim, p0*np.exp(-ylim/H)]
ylabs = ['height (km)', 'pressure (mb)']
yscales = ['height', 'pressure']
ylocators = [5, None]
x = np.linspace(*ylim, N)
xs = [x, 1000.0*np.exp(-x/H)]
y = np.cumsum((np.random.rand(len(x))-0.5), axis=0)
y = y - min(y)
colors = ['gray5', 'gray7']
ls = ['-', '--']
label = 'z = scale height = 7km, p = p$_{0}$/e = 368mb'
kw = dict(y=7, color='red', label=label, lw=2)
for i,ax in enumerate(axs):
    i = 1-i
    ax.plot(y, xs[i], color=colors[i], lw=2, ls=ls[i])
    ax.format(ylim=ylims[i], xlabel='quantity (units)', ylabel=ylabs[i],
              ylocator=ylocators[i], gridminor=True,
              suptitle='Profiles with pressure and height as the linear scale', abc=True)
    if i==0:
        h = ax.axhline(**kw)
    ax = ax.twinx()
    i = 1-i
    ax.format(ylim=ylims[i], ylabel=ylabs[i], yscale=yscales[i], ylocator=ylocators[i])
    if i==0:
        h = ax.axhline(**kw)
f.bottompanel.legend([h])
```

<!--
Configured ipython notebook.

```
{:.output_stream}
```
-->




<!--
{:.output_data_text}
```
(<matplotlib.axes._subplots.PanelAxesSubplot at 0x118e79748>,
 [<matplotlib.legend.Legend at 0x118fb98d0>])
```
-->



<!--
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_22_3.svg)


## Cutoff scales
Use so-called "cutoff scales" to create x/y axes with discrete cutoffs, or to
have x/y axes span different magnitudes across different parts of the axis.
Useful when you have data with large outliers.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
# plot.rc.fontname = 'Verdana'
f, axs = plot.figure(width=6, nrows=4, aspect=(5,1),
                     hspace=0.5,
                     sharey=False, sharex=False)
# Compression
ax = axs[0]
x = np.linspace(0,4*np.pi,1000)
xticks = plot.arange(0,12,1.0)
y = np.sin(x)
y2 = np.cos(x)
scales = [(3, np.pi), (0.3, 3*np.pi), (np.inf, np.pi, 2*np.pi), (5, np.pi, 2*np.pi)]
titles = ('Zoom out of left', 'Zoom into left', 'Discrete cutoff', 'Fast jump')
locators = [np.pi/3, np.pi/3, *([x*np.pi for x in plot.arange(0, 4, 0.25) if not (1 < x <= 2)] for i in range(2))]
for ax,scale,title,locator in zip(axs,scales,titles,locators):
    ax.plot(x, y, lw=3, color='blue7')
    ax.plot(x, y2, lw=3, color='red7')
    ax.format(xscale=('cutoff', *scale), title=title,
              xlim=(0,4*np.pi), ylabel='Wave amplitude', # note since 'spanning labels' turned on by default, only one label is drawn
              xformatter='pi', xlocator=locator,
              xtickminor=False, xgrid=True, ygrid=False)
```

<!--
Configured ipython notebook.
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_24_1.svg)


# Better colors
ProPlot provides several tools for creating plots with beautiful color palettes.

## New color names
This first plot shows newly registered colors from the [Open
Color](https://github.com/yeun/open-color) UI-design library. The second plot
shows the remaining registered colors, consisting of the standard ROYGBIV names,
"crayon" color names, and XKCD crowd-sourced color names. I limit the named
colors to those sufficiently distinct in the HCL colorspace (see below), to
eliminate redundant colors.


```python
import proplot as plot
plot.nbsetup()
f = plot.color_show(['open'])
```

<!--
Configured ipython notebook.
Saving to "/Users/ldavis/proplot/proplot/colors/colors_open.pdf".

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_27_1.svg)



```python
import proplot as plot
plot.nbsetup()
f = plot.color_show(nbreak=13)
```

<!--
Configured ipython notebook.
Saving to "/Users/ldavis/proplot/proplot/colors/colors_xkcd-crayons.pdf".

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_28_1.svg)


## Discrete colormaps
The below figure shows the newly regsistered discrete colormaps or "color
cycles" -- i.e., `ListedColormap`s, whose transitions are not meant to be
smooth. Any discrete colormap name can be used as the `cmap` argument in a
plotting command (e.g. `contourf`), and any smooth colormap name can be used as
the `cycler` argument in a plotting command (e.g. `plot`), or as the default
cycle `plot.rc.cycle`, using `cycle=('smooth_cmap', N)` where `N` indicates the
number of colors you wish to draw.


```python
import proplot as plot
plot.nbsetup()
f = plot.cycle_show()
```

<!--
Configured ipython notebook.
Saving to "/Users/ldavis/proplot/proplot/colors/cycles.pdf".

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_30_1.svg)


## Intro to colorspaces
My colormap generating tools, and some of the colormaps I provide by default,
are based on linear transitions for each channel in any of the following three
HSV-like colorspaces.

The **HCL colorspace** is a purely perceptually uniform colorspace, where colors
are broken down into "hue" (color, range 0-360), "chroma" (colorfulness, range
0-100), and "luminance" (brightness, range 0-100).

The problem is, many colors in the numeric range of this space are "imaginary"
(i.e. when converted to RGB, some channels exceed 1). We can "clip" the RGB
channels when this happens, or try a different approach: the HSLuv colorspace,
or the HPLuv colorspace.

The **HPLuv** colorspace scales 100 chroma to be the *minimum* max chroma across
*all hues for a given luminance*, and is hence more appropriate for multi-hue
colormaps. The **HSLuv** colorspace scales 100 chroma to be the *maximum
possible for a given hue and luminance*, and is hence more appropriate for
single-hue colormaps (crossing hues in this space make it more likely that bands
of higher absolute chroma are crossed; see the hue-luminance cross-section).

For more info, check out [this page](http://www.hsluv.org/comparison/).


```python
import proplot as plot
plot.nbsetup()
f = plot.colorspace_breakdown(luminance=50)
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_32_1.svg)



```python
import proplot as plot
plot.nbsetup()
f = plot.colorspace_breakdown(chroma=60)
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_33_1.svg)



```python
import proplot as plot
plot.nbsetup()
f = plot.colorspace_breakdown(hue=0)
```

<!--
Configured ipython notebook.
Adjusting gridspec.
Resetting rcparams.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_34_1.svg)



```python
import proplot as plot
plot.nbsetup()
plot.cmap_breakdown('NegPos')
```

<!--
Configured ipython notebook.

```
{:.output_stream}
```
-->

<!--
/Users/ldavis/anaconda3/lib/python3.6/site-packages/matplotlib/contour.py:1557: UserWarning: Warning: converting a masked element to nan.
  self.zmax = float(z.max())
/Users/ldavis/anaconda3/lib/python3.6/site-packages/matplotlib/contour.py:1558: UserWarning: Warning: converting a masked element to nan.
  self.zmin = float(z.min())

```
{:.output_stream}
```
-->

<!--
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_35_3.svg)



```python
import proplot as plot
plot.nbsetup()
plot.cmap_breakdown('Sunset')
```

<!--
Configured ipython notebook.

```
{:.output_stream}
```
-->

<!--
/Users/ldavis/anaconda3/lib/python3.6/site-packages/matplotlib/contour.py:1557: UserWarning: Warning: converting a masked element to nan.
  self.zmax = float(z.max())
/Users/ldavis/anaconda3/lib/python3.6/site-packages/matplotlib/contour.py:1558: UserWarning: Warning: converting a masked element to nan.
  self.zmin = float(z.min())

```
{:.output_stream}
```
-->

<!--
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_36_3.svg)


##  Smooth colormaps
By default, matplotlib comes packaged with every
[ColorBrewer2.0](http://colorbrewer2.org/) colormap. I've removed some outdated
"miscellaneous" colormaps that are packaged by default (see [this
reference](https://matplotlib.org/examples/color/colormaps_reference.html)), and
added the [cmOcean colormaps](https://matplotlib.org/cmocean/), and some pre-
packaged colormaps I generated with the `PerceptuallyUniformColormap` class,
which employs linear transitions for each channel in any of the perceptually
uniform colorpsaces. Note **every colormap can be referenced by its lower-case
name.**


```python
import proplot as plot
plot.nbsetup()
f = plot.cmap_show(31)
```

<!--
Configured ipython notebook.
Ignored colormaps: Wistia, afmhot, autumn, binary, bone, cool, coolwarm, copper, gist_gray, gist_heat, gist_yarg, gray, pink, seismic, spring, summer, winter, cividis, cubehelix, multi, bwr, coolwarm, seismic, CMRmap, brg, flag, gist_earth, gist_ncar, gist_rainbow, gist_stern, gnuplot, gnuplot2, hot, hsv, jet, nipy_spectral, ocean, prism, rainbow, terrain
Adjusting gridspec.
Saving to "/Users/ldavis/proplot/proplot/cmaps/colormaps.pdf".

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_38_1.svg)


## Cmap specification
This is one of the most versatile features offered by ProPlot. Colormaps can be
declared as gradations of a single color (e.g. `maroon`), registered colormap
names (e.g. `glacial`), registered color cycle names (e.g. `tropical`), lists of
colors (the `listmap` list below), and arbitrary linear transformations in HSL
space (the `flymap` dictionary below). Color cycles can be declared in the same
way as colormaps, optionally with an iterable `(cmap argument(s), N)` where `N`
indicates the number of colors you wish to draw.

To **concatenate** arbitrary colormaps, just pass an iterable containing the
previously described "colormap indicators" (e.g. `('C0', 'C1')` concatenates two
single-hue dark-light gradation colormaps). Note this means, if you want a
`ListedColormap` from a list of input colors, you must use `cmap=[['color1',
'color2', ...]]` insteaad of `cmap=['color1', 'color2', ...]`.

To **clip** colors in your colormap, use `left=n`, `right=m`, or `x=(n,m)` where
`n` and `m` are between 0 and 1.

To **build** a perceptually uniform colormap on-the-fly pass a dictionary with
keys `l[uminance]`, `c[hroma]`, and `h[ue]` (you can pass the first letters or
the full words). The corresponding values should indicate the range of hues,
luminances, and chromas across which you want your colormap to linearly vary.
Specify color name strings, and ProPlot will look up the corresponding channel
value for that color. Specify `'string+/-number'` to offset the channel value
for that color by any number. Use a `gamma>1` to make the colormap "linger" on
brighter/less colorful colors (i.e. the transitions will not be exactly linear).
Note **hues vary from 0 to 360**, while **luminance and chroma vary from 0 to
100**.


```python
import numpy as np
import proplot as plot
plot.nbsetup()
flymap = {'h':['blue-360','red'], 'l':[98, 20], 'space':'hpl', 'gamma':1.4}
listmap = ('light green', 'blue violet', 'sky blue', 'blue green', 'red violet')
cmaps  = ['maroon',     ('C0','C2'),    'tropical', 'glacial',         flymap,     [listmap], 'blood', 'blood']
cycles = [('maroon',N), ('C0','C2', N), 'tropical', ('glacial', N//2), (flymap, 5), listmap,  'blood', 'blood']
kws = [{}]*(len(cycles) - 1) + [{'left':0.3, 'right':0.9}] # clip colors on the last colormap
f, axs = plot.subplots(ncols=2, nrows=(len(cmaps)+1)//2,
                       axwidth=3, aspect=(5,4), share=3,
                       innerpanels_kw={'hspace':0.1, 'wwidth':0.8}, hspace=0.1,
                       innerpanels='r', innercolorbars='b')
# Lines
N = 12
lines = np.random.rand(20,N) - 0.5
lines = lines[:,:1] + lines.cumsum(axis=0) + np.arange(0,N)
ylim = (0,11)
scales = [0.1, 0.3, 0.5, 0.7]
for i,(ax,cmap,cycle,kw) in enumerate(zip(axs,cmaps,cycles,kws)):
    data = np.cos(np.sin(scales[i//2] * np.linspace(0,N,N)[None,:] * np.linspace(0,N,N)[:,None])) # psychadelic colors
    m = ax.contourf(data, cmap=cmap, cmap_kw=kw, levels=10)
    # ax.contour(data, colors='w', linewidths=0.5)
    ax.rightpanel.plot(lines, lw=2, cycle=cycle, cycle_kw=kw) # one for each line
    ax.rightpanel.format(ylocator='none', ylim=ylim)
    ax.bottompanel.colorbar(m, clocator='none')
axs.format(suptitle='Various ways to declare colormaps and cycles', abc=True, abcpos='il',
           xlim=None, xticks='none', ylim=ylim)
```

<!--
Configured ipython notebook.
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_40_1.svg)


## Other features
For any PerceptuallyUniformColormap, the chroma gamma (`gamma1`) and the
luminance gamma (`gamma2`) can be changed on-the-fly. For the former, large
numbers favor pale colors; for the latter, large numbers favor bright colors.
Thus this essentially gives the 'white' part of sequential/diverging cmaps more
emphasis.

Note I've also added support for pcolormesh *levels* and "extend" options (not
provided by default API). This is often very useful for interpreting physical
data with coarse resolution.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
f, axs = plot.subplots(ncols=3, nrows=2, innercolorbars='r',
                       hspace=0.3, wspace=0.2, aspect=1.2,
                       bspace=0.1)
data = np.random.rand(10,10).cumsum(axis=1)
def show(ax, cmap, gamma):
    m1 = ax.pcolormesh(data, cmap=cmap, cmap_kw={'gamma2':gamma}, levels=10, extend='both')
    ax.rightpanel.colorbar(m1, clocator='none')
    ax.format(title=f'gamma = {gamma}', xlabel='x axis', ylabel='y axis', suptitle='Varying gamma, and demo of new pcolor options')
cmap = 'verdant'
show(axs[0], cmap, 0.8)
show(axs[1], cmap, 1.0)
show(axs[2], cmap, 1.4)
cmap = 'fire'
show(axs[3], cmap, 0.8)
show(axs[4], cmap, 1.0)
show(axs[5], cmap, 1.4)
```

<!--
Configured ipython notebook.
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_42_1.svg)


I also enhanced the `plot` method to allow mapping colormap colors to each (x,y)
pair on a line. Use `interp=n` to interpolate `n` additional points between the
provided (x,y) pairs and corresponding colormap values, `values`.

Also demonstrated below is the "stacked colorbar", which is especially useful
where you have multiple colormaps on the same axes.


```python
import proplot as plot
import numpy as np
plot.nbsetup()
# Make a pretty spiral
N = 12
values = np.arange(1, N+1)
radii = np.linspace(1,0.2,N)
angles = np.linspace(0,4*np.pi,N)
# Figure
f, axs = plot.subplots(bottomcolorbar=True, ncols=2, wspace=0.35, aspect=1, axwidth=2.2, bwidth=0.8, span=False)
cmaps = [('blues', 'reds'), 'golden']
multipliers = [1.2, 1.4]
for i,(ax,cmap) in enumerate(zip(axs,cmaps)):
    x = radii*np.cos(multipliers[i]*angles)
    y = radii*np.sin(multipliers[i]*angles)
    m = ax.plot(x, y, cmap=cmap, values=values+i*12,
                linewidth=15, interp=1-i, cmap_kw={'left':i*0.05})
    ax.format(xlim=(-1,1), ylim=(-1,1), suptitle='Lines with smooth colormap gradations',
              xlabel='cosine angle', ylabel='sine angle')
    ax, c = f.bottompanel.colorbar(m,  space=0.37, i=i, n=2, locator=None, label=f'label {i}')
```

<!--
Configured ipython notebook.
Adjusting gridspec.

```
{:.output_stream}
```
-->


![svg](/tools/files/proplot_showcase_files/proplot_showcase_44_1.svg)

