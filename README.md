## Overview
This library provides extremely helpful plotting utilities to make `matplotlib` less painful. I recommend importing the utilities with `import pubplot as plot`. Most of the features below derive from the **`subplots`** command, a wrapper around the `pyplot` command of the same name.

Quick summary of features:

  * Shared and spanning axes: axis labels spanning subplots, shared axis tick labels.
  * Side panels: Outer panels spanning multiple rows/columns, suitable for shared legends or colorbars.
  * New and improved legend/colorbar: More flexible inputs, new features.
  * Sizing: Control over subplot axes ratios; set subplot spacing and panel sizes in inches (no more relative units).
  * Projection subplots: Integration with basemap/cartopy.
  * Colors: New colormaps, new colors, new color cycles, and functions for previewing them.

The **two most important** features are the `subplots` command and the `format` method. Use `subplots` to generate a scaffolding of axes and panels with rigid spacing, then use the new `format` method assigned to each axes instance to control the look of your axes in one line with a plethora of keyword arguments. Why do this?

  * To modify an axes property (e.g. an x-axis label) with the default API, you normally have to use a bunch of one-liner `pyplot` commands (or method calls on axes/axis objects). This can get repetitive and quite verbose, resulting in lost of copy-paste code.
  * Now, you can pass all of these settings to `format`. Instead of having to remember the name of the function, whether it's attached to `pyplot` or an object instance, and the order/names of the arguments, you just have to remember one thing -- the name of the keyword argument. The method also abstracts away some irritating inconsistencies -- now, There's Only One (obvious) Way To Do It.

### The `subplots` command
   * Generate grids of subplots:
     * Pass no arguments to generate default single-axes figure (1 row, 1 column).
     * Use `nrows` and `ncols` to generate simple subplot grids -- e.g. `nrows=2` creates two rows with one column.
     * To generate complex grids, pass a 2D array of numbers corresponding to unique subplots. Use zero to allot empty space. For example, `subplots(array=[[1,2],[1,3]])` creates a grid with one tall plot on the left,
     and two smaller plots on the right, while `subplots(array=[[1,1,1],[2,0,4]])` creates a grid with one long plot on top and two smaller plots on the bottom with a gap in the middle.
   * Returns two arguments: the figure, and a list of axes. If only one axes was drawn, returns the axes instead of a singleton list.
   * Precise control of layout: control figure width/height with kwargs `width/height`, vertical/horizontal space between axes with `wspace/hspace`, panel widths with `lwidth/cwidth/bwidth/rwidth`, inner-panel widths with `ihwidth/iwwidth`, inner-panel spacing with `ihspace/iwspace`, ratios for axes column widths (row heights) with `wratios/hratios`, ...
   * Automatic sizing: specify a width (height) with a subplot aspect ratio, and `subplot` will automatically determine the necessary height (width).
   * Integrated journal standards: for example, `width='ams1'` selects the smallest American Meteorological Society standard figure width.
### Useful features for multi-axes figures
   * Use `<bottom|right>panel=True` to allot space for panels spanning all columns (rows) on the bottom (right). Use `<bottom|right>panels=True` to allot space for one panel per column (row). Use `<bottom|right>panels=[n1,n2,...]` to allot space for panels that can span adjacent columns (rows) -- for example, if your subplot has 3 columns, passing `bottompanels=[1,2,2]` draws one panel for the first column and a second panel spanning the next two.
   * The above adds the `fig.<bottom|right>panel` attributes to figure object `fig`. These attributes are actually `SubplotSpec` instances. If your panel has more than one space, use `fig.<bottom|right>panel[n]` to access the nth space.
   * Use `fig.<bottom|right>panel.<legend|colorbar>` (or `fig.<bottom|right>panel[n].<legend|colorbar>`) to fill the `SubplotSpec`s with a legend or colorbar. Note: You must supply the legend command with a list of handles/supply the colorbar command with a mappable instance (e.g. `m = ax.contourf(...)`, `fig.bottompanel.colorbar(m, ...)`).
   * Outdated: Use `<bottom|right>colorbar=True` and  `fig.<bottom|right>colorbar.format(mappable, ...)` to create axes-spanning colorbars. Use `bottomlegend=True` and `fig.bottomlegend.format(handles, ...)` to create axes-spanning legends.
### The `format` command
   * To set axes properties, just pass kwargs to `format`. For example: `ax.format(xlabel="foo", ylabel="bar", title="foobar")`. Includes special kwargs that avoid using unpalatable underlying API.
   * Axis options:
   
      `xgrid`, `ygrid`,
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
      
    * Titling options:
    
      `suptitle`, `suptitlepos`, `title`, `titlepos`, `titlepad`, `titledict`,
      `abc`, `abcpos`, `abcformat`, `abcpad`, `abcdict`   
      
    * Mapping options:
    
      `oceans`, `coastlines`, `continents`,
      `latlabels`, `lonlabels`,
      `latlocator`, `latminorlocator`, `lonlocator`, `lonminorlocator` 
      
    * Axes canvas options:
    
       `hatch`, `color`
       
    * Example: draw major ticks every `2` units with `xlocator=2`, or specify a custom range with `xlocator=[0 2 4 8 16]`, instead of digging into `matplotlib.ticker.Locator` classes.
### Mapping toolkit integration
   * For projection subplots, specify `projection='<name>'` with either `package='basemap'` or `package='cartopy'`. Extra arguments to `subplot` will be passed to the `basemap.Basemap` and `cartopy.crs.Projection` classes (the relevant cartopy class will be selected based on the `'<name>'` string).
   * Control which subplots are projection subplots with `maps=[n1,n2,...]`, where numbers correspond to the subplot array number. Note that if axes numbers were not declared with `array`, the subplots are automatically numbered from 1 to n (row major).
   * `Basemap` instances are added as the attribute `m` their corresponding axes; create plots with (e.g.) `ax.m.contourf`. These instances are also overwritten to fix issues with `seams` on the edge of the map -- data will be circularly rolled and interpolated to map edges, so that seams are eliminated.
### Revised underlying issues with contour and pcolor commands
   * Flipped the unnatural default used by `pcolor` and `contour` functions: that `0`th dimension of the input array is `y`-axis, `1`st dimension is `x`-axis. More intuitive to enter array with `0`th dimension on `x`-axis.
   * The well-documented [white-lines-between-filled-contours](https://stackoverflow.com/q/8263769/4970632)nd [white-lines-between-pcolor-rectangles](https://stackoverflow.com/q/27092991/4970632) problems are fixed by automatically changing the edgecolors when `contourf`, `pcolor`, and `pcolormesh` are called.
### Enhanced settings management
   * Added the new `globals` command to change various settings. This command controls entries to the built-in `matplotlib.rcParams` dictionary along with the new, custom-built `matplotlib.rcExtras` dictionary.
   * Special arguments **integrate** settings across a wide range of pyplot objects -- for example, `globals(linewidth=1, small=8, large=9)` makes all axis back ground lines (axes spines, tick marks, colorbar edges, legend boxes, ...) have 1pt linewidth, makes the `small` size font (used for tick labels, axis legends, ...) 8pt, and makes the `large` size font (used for titles, 'abc' subplot labelling, ...) 9pt.
### Pretty colors
   * All colormaps from `.rgb` files in the `cmaps` folder can be called by name. For example, `hclBlue.rgb` can be used in a contour plot with `ax.contourf(x,y,z,cmap="hclBlue")`.
   * Added many new colors, including XKCD colors and Open Color web-design color palette. Added different color "cycles" (i.e. the automatic color order when drawing multiple lines), which can be set with `globals(cycle='<name>')`.
   * Use functions `cmapshow`, `colorshow`, and `cycleshow` to visualize available colormaps, named colors, and color palettes. Functions will automatically save PDFs in the package directory.
### Fonts
   * List of system fontnames available as the `fonts` variable under the imported module.
### Misc tools
   * New `arange` utility -- like `np.arange`, but **endpoint-inclusive**.
   * Default `Formatter` class for ticklabels renders numbers into the style you'll want 90% of the time. Also use `LatFormatter` or `LonFormatter` for coordinate axes. Former can be used to draw sine-weighted latitude axes.
   * New `Normalize` class allows colormap data limits to have `cmin != cmax`. For example, you can create contours with levels `arange(-20,50,10)` and keep the midpoint of the colormap at `0`, with either ends fully saturated.

