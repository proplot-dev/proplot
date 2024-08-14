import proplot as plt, numpy as np
import pytest


@pytest.mark.mpl_image_compare
def test_geographic_single_projection():
    fig = plt.figure(refwidth=3)
    axs = fig.subplots(nrows=2, proj="robin", proj_kw={"lon_0": 180})
    # proj = pplt.Proj('robin', lon_0=180)
    # axs = pplt.subplots(nrows=2, proj=proj)  # equivalent to above
    axs.format(
        suptitle="Figure with single projection",
        land=True,
        latlines=30,
        lonlines=60,
    )
    return fig


@pytest.mark.mpl_image_compare
def test_geographic_multiple_projections():
    fig = plt.figure()
    # Add projections
    gs = plt.GridSpec(ncols=2, nrows=3, hratios=(1, 1, 1.4))
    for i, proj in enumerate(("cyl", "hammer", "npstere")):
        ax1 = fig.subplot(gs[i, 0], proj=proj, basemap=True)  # basemap
        ax2 = fig.subplot(gs[i, 1], proj=proj)  # cartopy

    # Format projections
    fig.format(
        land=True,
        suptitle="Figure with several projections",
        toplabels=("Basemap projections", "Cartopy projections"),
        toplabelweight="normal",
        latlines=30,
        lonlines=60,
        lonlabels="b",
        latlabels="r",  # or lonlabels=True, labels=True, etc.
    )
    fig.subplotgrid[-2:].format(
        latlines=20, lonlines=30
    )  # dense gridlines for polar plots
    plt.rc.reset()
    return fig


@pytest.mark.mpl_image_compare
def test_drawing_in_projection_without_globe():
    # Fake data with unusual longitude seam location and without coverage over poles
    offset = -40
    lon = plt.arange(offset, 360 + offset - 1, 60)
    lat = plt.arange(-60, 60 + 1, 30)
    state = np.random.RandomState(51423)
    data = state.rand(len(lat), len(lon))

    globe = False
    string = "with" if globe else "without"
    gs = plt.GridSpec(nrows=2, ncols=2)
    fig = plt.figure(refwidth=2.5)
    for i, ss in enumerate(gs):
        ax = fig.subplot(ss, proj="kav7", basemap=(i % 2))
        cmap = ("sunset", "sunrise")[i % 2]
        if i > 1:
            ax.pcolor(lon, lat, data, cmap=cmap, globe=globe, extend="both")
        else:
            m = ax.contourf(lon, lat, data, cmap=cmap, globe=globe, extend="both")
            fig.colorbar(m, loc="b", span=i + 1, label="values", extendsize="1.7em")
    fig.format(
        suptitle=f"Geophysical data {string} global coverage",
        toplabels=("Cartopy example", "Basemap example"),
        leftlabels=("Filled contours", "Grid boxes"),
        toplabelweight="normal",
        leftlabelweight="normal",
        coast=True,
        lonlines=90,
        abc="A.",
        abcloc="ul",
        abcborder=False,
    )
    return fig


@pytest.mark.mpl_image_compare
def test_drawing_in_projection_with_globe():
    # Fake data with unusual longitude seam location and without coverage over poles
    offset = -40
    lon = plt.arange(offset, 360 + offset - 1, 60)
    lat = plt.arange(-60, 60 + 1, 30)
    state = np.random.RandomState(51423)
    data = state.rand(len(lat), len(lon))

    globe = True
    string = "with" if globe else "without"
    gs = plt.GridSpec(nrows=2, ncols=2)
    fig = plt.figure(refwidth=2.5)
    for i, ss in enumerate(gs):
        ax = fig.subplot(ss, proj="kav7", basemap=(i % 2))
        cmap = ("sunset", "sunrise")[i % 2]
        if i > 1:
            ax.pcolor(lon, lat, data, cmap=cmap, globe=globe, extend="both")
        else:
            m = ax.contourf(lon, lat, data, cmap=cmap, globe=globe, extend="both")
            fig.colorbar(m, loc="b", span=i + 1, label="values", extendsize="1.7em")
    fig.format(
        suptitle=f"Geophysical data {string} global coverage",
        toplabels=("Cartopy example", "Basemap example"),
        leftlabels=("Filled contours", "Grid boxes"),
        toplabelweight="normal",
        leftlabelweight="normal",
        coast=True,
        lonlines=90,
        abc="A.",
        abcloc="ul",
        abcborder=False,
    )
    return fig
