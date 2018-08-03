## Overview
This library provides extremely helpful plotting utilities to make `matplotlib` less painful,
and provides sort-of-helpful-but-still-working-on-it numerical utilities for analysing data.

## Plotting Utilities
Everything here is contained in `plots`. Usually I import these as `from pyfuncs import plot`. Most of the features below derive from the **`subplots`** command, a wrapper around the `pyplot` version.

Quick summary of features:

  * Total control over subplot aspect ratios, figure sizes in inches, spacing between subplots and around edges in inches.
  * Tools to generate panels around `edges` of primary subplots, suitable for global legends/global colorbars/
  * Axis labels spanning subplots and no more duplicate tick labels.
  * Inner subplot panels, e.g. for displaying averages across the x/y axis.
  * Integration with basemap/cartopy.

The **most notable** feature is the `format` method added to axes instances. Now, you control the look of your axes with a single function, that accepts a plethora of kwarg pairs. Why do this?

  * With the default API, to modify an axes property, you generally have to remember **three things**: 1) the object containing the method you want to use (e.g. axes vs. an individual axis), 2) the name of the method, and 3) the keyword-arguments and their meaning. It can get quite messy, quite verbose, and inevitably results in lots of copy-paste code!
  * Now, you just have to remember **one** thing: the name of the keyword-argument passed to `ax.format()`.

For more information on the `plot` features, see below.

### Improved functionality of `subplots` command
   * Much easier generation of figures with multiple axes.
     * Pass no arguments to generate default single-axes figure (1 row, 1 column).
     * Pass the `nrows` and/or `ncols` kwarg to generate simple grids of axes -- e.g. `nrows=2` creates two rows with one column.
     * Pass an array with numbers corresponding to unique subplots to generate complex grids -- e.g. `py.subplots([[1,2],[1,3]])` creates a grid with one tall plot on the left,
     and two smaller plots on the right, while `py.subplots([[1,1,1],[2,3,4]])` creates a grid with one long plot on top and three smaller plots on the bottom.
   * Much more precise sizing control -- control in **inches** the figure width and/or height, vertical/horizontal space between axes, ratios of axes row heights/axes column widths, extra spacing around outside of leftmost/rightmost columns and bottom/top rows. Allows making publication-quality graphics without the publishers needing to re-scale them.
   * Added many new features -- see below.
### Much-wanted features for multi-axes figures
   * Use `bottomcolorbar=True` and `rightcolorbar=True` to create special axes-spanning colorbars on your plots. Can be accessed as member of `Figure` instance.
   * Use `bottomlegend=True` to create axes-spanning legend at the bottom. Can be accessed as member of `Figure` instance.
   * Exact widths and spacings can be controlled in inches.
### Integrated mapping toolkits seamlessly
   * Can pass `projection=<name>` with either `package='basemap'` or `package='cartopy'`. Pass extra map arguments to the `subplots` command directly.
   * This creates axes grids, with each axes a map. Power to choose between cartopy and basemap.
### Added new `format` method-style commands to axes instances generated with `subplots`
   * To customize an axes, most important properties can be passed as kwargs. For example, format an axes with `ax.format(xlabel="foo", ylabel="bar", title="foobar")`.
   * Included some special kwargs that avoid using unpalatable underlying API. For example, draw major ticks every `2` units with `xlocator=2`, or specify a custom range with `xlocator=[0 2 4 8 16]`.
   * Special `format` functions created for `bottomlegend` and `bottomcolorbar`/`rightcolorbar`.
### Revised some underlying issues by overwriting existing methods
   * Flipped the unnatural default `0`th dimension of array is `y`-axis, `1`st dimension of array is `x`-axis. More intuitive to enter array with `0`th dimension on `x`-axis.
   * The well-documented [white-lines-between-filled-contours](https://stackoverflow.com/q/8263769/4970632)nd [white-lines-between-pcolor-rectangles](https://stackoverflow.com/q/27092991/4970632) problems are fixed by automatically changing the edgecolors when `contourf`, `pcolor`, and `pcolormesh` are called.
   * Wrappers around the special `Basemap` plotting utilities are added as methods to axes. Can call them with `mcontourf`, `mpcolor`, `mpcolormesh`, and `mcontour`.
   These roll data to appropriate longitude range and fix seams on the map edges.
### More flexible settings management compared to `mpl.rcParams`
   * Can change settings for various axes objects with `py.setup(dict1, dict2)`; for example `py.setup({'ticklabels':size=5})`.
   * Settings are applied whenever the `format` method is called, kept constant through entire workspace.
   * Settings can be accessed at any time under `py.settings.ticklabels`, `py.settings.title`, etc.
### Font management
   * List of system fontnames available as the `fonts` variable under the imported module.
### Colormap management
   * On import, colormaps are automatically added from any `.rgb` files in the `cmaps` folder, and can be called by name. For example, `hclBlue.rgb` can be used in a contour plot with `ax.contourf(x,y,z,cmap="hclBlue")`.
      * Files downloaded from the HCL Wizard should be prefixed with `hcl`.
      * Other formats must be coded into the `plots.py` function to be read properly
   * Created function `cmapshow` to show off current colormaps in nice plot, automatically saved as `colormaps.pdf` in repository.
### Extra tools
   * Created new `arange` utility. Like `np.arange`, but this is always **endpoint-inclusive**.
   * Created `Formatter` class for ticklabels that renders numbers into the style you'll want 90% of the time.
   * Created `Normalize` class for colored content that allows for the colormap data limits to have `cmin != cmax`. For example, you can create contours with levels `np.arange(-20,51,10)` and keep the midpoint of the colormap at `0`, with either ends fully saturated.

## Numerical Utilities
Need to document these better, decide if this duplicates anything.


