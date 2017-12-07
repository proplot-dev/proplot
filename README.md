## Overview
This library provides extremely helpful plotting utilities to make `matplotlib` less painful,
and sort-of-helpful-but-still-working-on-it numerical utilities for analysing data.

## Plotting Utilities

Everything here is contained in `plots.py`. Usually I import these as `import pyfuncs.plots as py`. Current features:
### Improved `subplots` command
   * Instead of `plt.subplot`, use `py.subplot`.
   * Can pass `nrows`/`ncols` keyword arguments for simple grids.
   * Can pass an array with numbers corresponding to unique subplots for complex grids -- e.g. `py.subplots([[1,2],[1,3]])` creates a grid with one tall plot on the left,
   and two smaller plots on the right.
   * Control spaces between plots and height/width ratios with `wratios`, `hratios`, `wspace`, `hspace`.
### Much-wanted features for multi-axes figures
   * Use `bottomcolorbar=True` and `rightcolorbar=True` to create special axes-spanning colorbars on your plots. Can be accessed as member of `Figure` instance.
   * Use `bottomlegend=True` to create axes-spanning legend at the bottom. Can be accessed as member of `Figure` instance.
### Integrated mapping toolkits seamlessly
   * Can pass `projection=<name>` with either `package='basemap'` or `package='cartopy'`. Pass extra map arguments to the `subplots` command directly.
   * This creates axes grids, with each axes a map. Power to choose between cartopy and basemap.
### Added new `format` method-style commands to axes instances generated with `subplots`
   * To customize an axes, most important properties can be passed as kwargs. For example, format an axes with `ax.format(xlabel="foo", ylabel="bar", title="foobar")`.
   * Included some special kwargs that avoid using unpalatable underlying API. For example, draw major ticks every `2` units with `xlocator=2`, or specify a custom range with `xlocator=[0 2 4 8 16]`.
   * Special `format` functions created for `bottomlegend` and `bottomcolorbar/`rightcolorbar`.
### Revised some underlying issues by overwriting existing methods
   * Flipped the unnatural default `0`th dimension of array is `y`-axis, `1`st dimension of array is `x`-axis. More intuitive to enter array with `0`th dimension on `x`-axis.
   * The well-documented white-lines-between-my-solid-contours problem is fixes by automatically changing the edgecolor when `contourf` is called.
   * Wrappers around the special `Basemap` plotting utilities are added as methods to axes. Can call them with `mcontourf`, `mpcolor`, `mpcolormesh`, and `mcontour`.
   These roll data to appropriate longitude range and fix seams on the map edges.
### Extra tools
   * Created new `arange` utility. Like `np.arange`, but this is always **endpoint-inclusive**.
   * Created `Formatter` class for ticklabels that renders numbers into the style you'll want 90% of the time.
   * Created `Normalize` class for colored content that allows for the colormap data limits to have `cmin != cmax`. For example, you can create contours with levels `np.arange(-20,51,10)` and keep the midpoint of the colormap at `0`, with either ends fully saturated.
### Better colormap management
   * On import, colormaps are automatically added from any `.rgb` files in the `cmaps` folder, and can be called by name. For example, `hclBlue.rgb` can be used in a contour plot with `ax.contourf(x,y,z,cmap="hclBlue")`.
      * Files downloaded from the HCL Wizard should be prefixed with `hcl`.
      * Other formats must be coded into the `plots.py` function to be read properly
   * Created function `cmapshow` to show off current colormaps in nice plot, automatically saved as `colormaps.pdf` in repository.

## Numerical Utilities
Need to document these better, decide if this duplicates anything.
