#!/usr/bin/env python3
# import proplot as pplt
import numpy as np, pandas as pd, proplot as pplt
import pytest


@pytest.mark.mpl_image_compare
def test_statistical_boxplot():
    # Sample data
    N = 500
    state = np.random.RandomState(51423)
    data1 = state.normal(size=(N, 5)) + 2 * (state.rand(N, 5) - 0.5) * np.arange(5)
    data1 = pd.DataFrame(data1, columns=pd.Index(list("abcde"), name="label"))
    data2 = state.rand(100, 7)
    data2 = pd.DataFrame(data2, columns=pd.Index(list("abcdefg"), name="label"))

    # Figure
    fig, axs = pplt.subplots([[1, 1, 2, 2], [0, 3, 3, 0]], span=False)
    axs.format(abc="A.", titleloc="l", grid=False, suptitle="Boxes and violins demo")

    # Box plots
    ax = axs[0]
    obj1 = ax.box(data1, means=True, marker="x", meancolor="r", fillcolor="gray4")
    ax.format(title="Box plots")

    # Violin plots
    ax = axs[1]
    obj2 = ax.violin(data1, fillcolor="gray6", means=True, points=100)
    ax.format(title="Violin plots")

    # Boxes with different colors
    ax = axs[2]
    ax.boxh(data2, cycle="pastel2")
    ax.format(title="Multiple colors", ymargin=0.15)
    return fig


@pytest.mark.mpl_image_compare
def test_panel_dist():
    # Sample data
    N = 500
    state = np.random.RandomState(51423)
    x = state.normal(size=(N,))
    y = state.normal(size=(N,))
    bins = pplt.arange(-3, 3, 0.25)

    # Histogram with marginal distributions
    fig, axs = pplt.subplots(ncols=2, refwidth=2.3)
    axs.format(
        abc="A.",
        abcloc="l",
        titleabove=True,
        ylabel="y axis",
        suptitle="Histograms with marginal distributions",
    )
    colors = ("indigo9", "red9")
    titles = ("Group 1", "Group 2")
    for ax, which, color, title in zip(axs, "lr", colors, titles):
        ax.hist2d(
            x,
            y,
            bins,
            vmin=0,
            vmax=10,
            levels=50,
            cmap=color,
            colorbar="b",
            colorbar_kw={"label": "count"},
        )
        color = pplt.scale_luminance(color, 1.5)  # histogram colors
        px = ax.panel(which, space=0)
        px.histh(y, bins, color=color, fill=True, ec="k")
        px.format(grid=False, xlocator=[], xreverse=(which == "l"))
        px = ax.panel("t", space=0)
        px.hist(x, bins, color=color, fill=True, ec="k")
        px.format(grid=False, ylocator=[], title=title, titleloc="l")
    return fig
