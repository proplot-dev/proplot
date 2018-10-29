## Overview
This library provides helpful and versatile plotting utilities to make the process of crafting publication-quality graphics with `matplotlib` smoother. I recommend importing the package with
```
import pubplot as plot
```
Most of the features below derive from the **`subplots`** command, a wrapper around the `pyplot` command of the same name.

Quick summary of features:

  * Shared and spanning axes: axis labels spanning subplots, shared axis tick labels.
  * Side panels: Outer panels spanning multiple rows/columns, suitable for shared legends or colorbars.
  * New and improved legend/colorbar: More flexible inputs, new features.
  * Sizing: Control over subplot axes ratios; set subplot spacing and panel sizes in inches (no more relative units).
  * Projection subplots: Integration with basemap/cartopy.
  * Colors: New colormaps, new colors, new color cycles, and functions for previewing them.

The **two most important** features are the `subplots` command and the `format` method. Use `subplots` to generate a scaffolding of axes and panels with rigid spacing, then use the new `format` method assigned to each axes instance to control the look of your axes in one line with a plethora of keyword arguments. Why do this?

  * To modify an axes property (e.g. an x-axis label) with the default API, you normally have to use a bunch of one-liner `pyplot` commands (or method calls on axes/axis objects). This can get repetitive and quite verbose, resulting in lots of copy-paste code.
  * Now, you can pass all of these settings to `format`. Instead of having to remember the name of the function, whether it's attached to `pyplot` or an object instance, and the order/names of the arguments, you just have to remember one thing -- the name of the keyword argument. The method also abstracts away some irritating inconsistencies and redundancies -- now, There's Only One (obvious) Way To Do It.

## Installation
```
pip install git+https://github.com/lukelbd/pubplot.git#egg=pubplot
```

## Documentation
### The `subplots` command
   * Generate grids of subplots:
     * Pass no arguments to generate default single-axes figure (1 row, 1 column).
     * Use `nrows` and `ncols` to generate simple subplot grids -- e.g. `nrows=2` creates two rows with one column.
     * To generate complex grids, pass a 2D array of numbers corresponding to unique subplots. Use zero to allot empty space. For example, `subplots(array=[[1,2],[1,3]])` creates a grid with one tall plot on the left,
     and two smaller plots on the right, while `subplots(array=[[1,1,1],[2,0,4]])` creates a grid with one long plot on top and two smaller plots on the bottom with a gap in the middle.
   * Returns two arguments: the figure, and a list of axes. If only one axes was drawn, returns the axes instead of a singleton list.
   * Precise control of layout: control figure width/height with kwargs `width/height`, vertical/horizontal space between axes with `wspace/hspace`, panel widths with `bwidth/rwidth`, inner-panel widths with `ihwidth/iwwidth`, inner-panel spacing with `ihspace/iwspace`, ratios for axes column widths (row heights) with `wratios/hratios`, ...
   * Automatic sizing: specify a width (height) with a subplot *aspect ratio* (e.g. `aspect=1`, the default), and `subplot` will automatically determine the necessary height (width) required to preserve that aspect ratio.
   * Integrated journal standards: for example, `width='ams1'` selects the smallest American Meteorological Society standard figure width.
### Inner and outer "panels"
Use `[bottom|right]panel=True` to allot space for panels spanning all columns (rows) on the bottom (right). Use `[bottom|right]panels=True` to allot space for one panel per column (row). Use `[bottom|right]panels=[n1,n2,...]` to allot space for panels that can span adjacent columns (rows) -- for example, if your subplot has 3 columns, passing `bottompanels=[1,2,2]` draws one panel for the first column and a second panel spanning the next two.

The above adds the `fig.[bottom|right]panel` attributes to figure object `fig`. These attributes are actually `SubplotSpec` instances. If your panel has more than one space, use `fig.[bottom|right]panel[n]` to access the nth space.

Convenience feature: `[bottom|right][colorbar|legend][s]=True` to modify the panel width and spacing to be *suitable for colorbars/legends*, e.g. `rightcolorbar=True`.

   * Use `fig.[bottom|right]panel.[legend|colorbar]` to fill the `SubplotSpec` with a legend or colorbar. Note: You must supply the legend command with a list of handles/supply the colorbar command with a mappable instance (e.g. `m = ax.contourf(...)`, `fig.bottompanel.colorbar(m, ...)`).
   * Use `fig.[bottom|right]panel` with any other plotting method (`plot`, `contourf`, etc.) to fill the `SubplotSpec` with an axes, then draw stuff on those axes.
### The `format` command
   * To set axes properties, just pass kwargs to `format`. For example: `ax.format(xlabel="foo", ylabel="bar", title="foobar")`. Includes special kwargs that avoid using unpalatable underlying API.
   * Axis options:  `xgrid`, `ygrid`,
      `xdates`, `ydates`,
      `xtickminor`, `ytickminor`, `xgridminor`, `ygridminor`,
      `xspineloc`, `yspineloc`,
      `xtickloc`, `ytickloc`,
      `xtickdir`, `ytickdir`,
      `xticklabeldir`, `yticklabeldir`,
      `xtickrange`, `ytickrange`,
      `xlim`, `ylim`, `xscale`, `yscale`, `xscale_kwargs`, `yscale_kwargs`,
      `xreverse`, `yreverse`,
      `xlabel`, `ylabel`,
      `xlocator`, `xminorlocator`, `ylocator`, `yminorlocator`,
      `xformatter`, `yformatter`    
   * Titling options: `suptitle`, `suptitlepos`, `title`, `titlepos`, `titlepad`, `titledict`,
      `abc`, `abcpos`, `abcformat`, `abcpad`, `abcdict`        
   * Mapping options: `oceans`, `coastlines`, `continents`,
      `latlabels`, `lonlabels`,
      `latlocator`, `latminorlocator`, `lonlocator`, `lonminorlocator` 
   * Axes canvas options: `hatch`, `color`      
   * Example: draw major ticks every `2` units with `xlocator=2`, or specify a custom range with `xlocator=[0 2 4 8 16]`, instead of digging into `matplotlib.ticker.Locator` classes.
### Mapping toolkit integration
   * For projection subplots, specify `projection='name'` with either `package='basemap'` or `package='cartopy'`. Extra arguments to `subplot` will be passed to the `basemap.Basemap` and `cartopy.crs.Projection` classes (the relevant cartopy class will be selected based on the `'name'` string).
   * Control which subplots are projection subplots with `maps=[n1,n2,...]`, where numbers correspond to the subplot array number. Note that if axes numbers were not declared with `array`, the subplots are automatically numbered from 1 to n (row major).
   * Access basemap plotting utilities directly as an axes method, thanks to the `BasemapAxes` subclass. Several plotting methods are also overridden to fix issues with the "seam" on the edge of the map (data is circularly rolled and interpolated to map edges).
### New x/y axis scales, tick formatters, and tick locators
   * Added scale for **sine-weighted** and **inverse-weighted** x or y-axes. Invoke with `[x|y]scale='sine'` and `[x|y]scale='inverse'`. The former is useful for plots against geographic latitude, the latter is useful where you wish to have both **wavenumber and wavelength** labeled on the opposite spines.
   * Added arbitrary scale factory that can create scales with custom cutoffs.
   * The new default `Formatter` class for ticklabels renders numbers into the style you'll want 90% of the time.
   * Pass `formatter='[lat|deglat|lon|deglon|deg]'` to format axis labels with cardinal direction indicators or degree symbols (as denoted by the names).
   * Pass `locator='[string]'` to use any of the `matplotlib.ticker` locators, e.g. `locator='month'` or `locator='log'`.
### Revised underlying issues with contour and pcolor commands
   * Flipped the unnatural default used by `pcolor` and `contour` functions: that `0`th dimension of the input array is `y`-axis, `1`st dimension is `x`-axis. More intuitive to enter array with `0`th dimension on `x`-axis.
   * The well-documented [white-lines-between-filled-contours](https://stackoverflow.com/q/8263769/4970632)nd [white-lines-between-pcolor-rectangles](https://stackoverflow.com/q/27092991/4970632) problems are fixed by automatically changing the edgecolors when `contourf`, `pcolor`, and `pcolormesh` are called.
### Enhanced settings management
   * Added the new `rc_configurator` class suitable for changing global settings. An instance named `rc` is created when `pubplot` is imported, and it can be used to change built-in `matplotlib.rcParams` settings, a few custom "`rcSpecial`" settings, and some special "global" settings that modify several other settings at once.
   * Example: Use `plot.rc['linewidth'] = 2` or `plot.rc.linewidth = 2` to increase the thickness of axes spines, major tick marks, and minor tick marks. Use `plot.rc['color'] = 'red'` or `plot.rc.color = 'red'` to make all spines, tick marks, tick labels, and axes labels red. There is also `plot.rc['small']` and `plot.rc['large']` to control axes font sizes. Update any arbitrary `rcParams` setting with e.g. `plot.rc['legend.frameon'] = False`.
   * Use `plot.rc.reset()` to reset everything to the initial state. This is also called every time a figure's `draw` method is invoked (e.g. when a figure is rendered by the matplotlib backend or saved to file).
   
### Colormaps, color cycles, and color names
   * Added **new colormap class** analagous to `LinearSegmentedColormap`, called `PerceptuallyUniformColormap`, which interpolates through hue, saturation, and luminance space (with hues allowed to vary circularly). Choose from either of 4 HSV-like colorspaces: classic HSV, perceptually uniform HCL, or HSLuv/HPLuv (which are forms of HCL adapted for convenient use with colormaps; see [this link](http://www.hsluv.org/comparison/)).
   * Generate perceptually uniformly varying colormaps on-the-fly by passing a **dictionary** to any plotting function that accepts the `cmap`keyword argument -- for example, `ax.contourf(..., cmap=dict(h=hues, l=lums, s=sats))`. The arguments can be lists of numbers or **color strings**, in which case the corresponding channel value (hue, saturation, or luminance) for that color will be looked up and applied.
   * Create single-hue colormaps on-the-fly by passing a string that looks like `cmap='name_[mode]'`, where `name` is any registered color string (the corresponding hue will be looked up). The colormap mode can be either of `l` (vary color from its current lightness to near-white), `d` (vary color/chroma from its current lightness/saturation to pure gray), `w` (as in `l`, but use pure white), and `b` (as in `d`, but use pure black). 
   * The `.rgb` and `.hex` files in the `cmaps` folder will be loaded and registered as `LinearSegmentedColormaps`. 
   , `Sunset.rgb` can be used in a contour plot with `ax.contourf(x,y,z,cmap='sunset')` (note cmap selection is also case insensitive now).
   * Added several new color "cycles" (i.e. the automatic color order when drawing multiple lines). Cycler can be set with `plot.rc.cycle = 'name'` or by passing `cycle='name'` to any command that plots lines/patches (`plot`, `bar`, etc.).
   * The definition of "colormaps" and "color cyclers" is now fluid -- all color cycles are defined as `ListedColormaps`, and cycles can be generated on the fly from `LinearSegmentedColormaps` by just specifying the colormap name as the cycler. 
   * Registered many new colors names, including XKCD colors, crayon colors, and Open Color web-design color palette.
   * Use functions `cmapshow`, `colorshow`, and `cycleshow` to visualize available colormaps, named colors, and color cycles. Functions will automatically save PDFs in the package directory.
### Miscellaneous
   * New `plot.arange` utility -- like `np.arange`, but **endpoint-inclusive**. This is useful for, e.g., designating tick positions and contour levels.
