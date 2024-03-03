import pytest

import proplot as plt, numpy as np
from matplotlib.testing import setup


@pytest.fixture()
def setup_mpl():
    setup()
    plt.clf()


@pytest.mark.mpl_image_compare
def test_standardized_input():
    # Sample data
    state = np.random.RandomState(51423)
    x = y = np.array([-10, -5, 0, 5, 10])
    xedges = plt.edges(x)
    yedges = plt.edges(y)
    data = state.rand(y.size, x.size)  # "center" coordinates
    lim = (np.min(xedges), np.max(xedges))

    with plt.rc.context({"cmap": "Grays", "cmap.levels": 21}):
        # Figure
        fig = plt.figure(refwidth=2.3, share=False)
        axs = fig.subplots(ncols=2, nrows=2)
        axs.format(
            xlabel="xlabel",
            ylabel="ylabel",
            xlim=lim,
            ylim=lim,
            xlocator=5,
            ylocator=5,
            suptitle="Standardized input demonstration",
            toplabels=("Coordinate centers", "Coordinate edges"),
        )

        # Plot using both centers and edges as coordinates
        axs[0].pcolormesh(x, y, data)
        axs[1].pcolormesh(xedges, yedges, data)
        axs[2].contourf(x, y, data)
        axs[3].contourf(xedges, yedges, data)
    fig.show()
    return fig


@pytest.mark.mpl_image_compare
def test_inbounds_data():
    # Sample data
    cmap = "turku_r"
    state = np.random.RandomState(51423)
    N = 80
    x = y = np.arange(N + 1)
    data = 10 + (state.normal(0, 3, size=(N, N))).cumsum(axis=0).cumsum(axis=1)
    xlim = ylim = (0, 25)

    # Plot the data
    fig, axs = plt.subplots(
        [[0, 1, 1, 0], [2, 2, 3, 3]],
        wratios=(1.3, 1, 1, 1.3),
        span=False,
        refwidth=2.2,
    )
    axs[0].fill_between(
        xlim,
        *ylim,
        zorder=3,
        edgecolor="red",
        facecolor=plt.set_alpha("red", 0.2),
    )
    for i, ax in enumerate(axs):
        inbounds = i == 1
        title = f"Restricted lims inbounds={inbounds}"
        title += " (default)" if inbounds else ""
        ax.format(
            xlim=(None if i == 0 else xlim),
            ylim=(None if i == 0 else ylim),
            title=("Default axis limits" if i == 0 else title),
        )
        ax.pcolor(x, y, data, cmap=cmap, inbounds=inbounds)
    fig.format(
        xlabel="xlabel",
        ylabel="ylabel",
        suptitle="Default vmin/vmax restricted to in-bounds data",
    )
    fig.show()
    return fig


@pytest.mark.mpl_image_compare
def test_colorbar():
    # Sample data
    state = np.random.RandomState(51423)
    data = 10 + state.normal(0, 1, size=(33, 33)).cumsum(axis=0).cumsum(axis=1)

    # Figure
    fig, axs = plt.subplots([[1, 1, 2, 2], [0, 3, 3, 0]], ref=3, refwidth=2.3)
    axs.format(yformatter="none", suptitle="Discrete vs. smooth colormap levels")

    # Pcolor
    axs[0].pcolor(data, cmap="viridis", colorbar="l")
    axs[0].set_title("Pcolor plot\ndiscrete=True (default)")
    axs[1].pcolor(data, discrete=False, cmap="viridis", colorbar="r")
    axs[1].set_title("Pcolor plot\ndiscrete=False")

    # Imshow
    m = axs[2].imshow(data, cmap="oslo", colorbar="b")
    axs[2].format(title="Imshow plot\ndiscrete=False (default)", yformatter="auto")
    fig.show()
    return fig
