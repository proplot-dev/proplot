#!/usr/bin/env python3
"""
Test 1D plotting overrides.
"""
import numpy as np
import numpy.ma as ma
import pandas as pd

import proplot as pplt
import pytest

state = np.random.RandomState(51423)


@pytest.mark.mpl_image_compare
def test_auto_reverse():
    """
    Test enabled and disabled auto reverse.
    """
    x = np.arange(10)[::-1]
    y = np.arange(10)
    z = state.rand(10, 10)
    fig, axs = pplt.subplots(ncols=2, nrows=3, share=0)
    # axs[0].format(xreverse=False)  # should fail
    axs[0].plot(x, y)
    axs[1].format(xlim=(0, 9))  # manual override
    axs[1].plot(x, y)
    axs[2].plotx(x, y)
    axs[3].format(ylim=(0, 9))  # manual override
    axs[3].plotx(x, y)
    axs[4].pcolor(x, y[::-1], z)
    axs[5].format(xlim=(0, 9), ylim=(0, 9))  # manual override!
    axs[5].pcolor(x, y[::-1], z)
    fig.format(suptitle="Auto-reverse test", collabels=["reverse", "fixed"])
    return fig


@pytest.mark.mpl_image_compare
def test_cmap_cycles():
    """
    Test sampling of multiple continuous colormaps.
    """
    cycle = pplt.Cycle(
        "Boreal",
        "Grays",
        "Fire",
        "Glacial",
        "yellow",
        left=[0.4] * 5,
        right=[0.6] * 5,
        samples=[3, 4, 5, 2, 1],
    )
    fig, ax = pplt.subplots()
    data = state.rand(10, len(cycle)).cumsum(axis=1)
    data = pd.DataFrame(data, columns=list("abcdefghijklmno"))
    ax.plot(data, cycle=cycle, linewidth=2, legend="b")
    return fig


@pytest.mark.mpl_image_compare
def test_column_iteration():
    """
    Test scatter column iteration.
    """
    fig, axs = pplt.subplots(ncols=2)
    axs[0].plot(state.rand(5, 5), state.rand(5, 5), lw=5)
    axs[1].scatter(
        state.rand(5, 5), state.rand(5, 5), state.rand(5, 5), state.rand(5, 5)
    )
    return fig


@pytest.mark.skip("TODO")
@pytest.mark.mpl_image_compare
def test_bar_stack():
    """
    Test bar and area stacking.
    """
    # TODO: Add test here


@pytest.mark.mpl_image_compare
def test_bar_width():
    """
    Test relative and absolute widths.
    """
    fig, axs = pplt.subplots(ncols=3)
    x = np.arange(10)
    y = state.rand(10, 2)
    for i, ax in enumerate(axs):
        ax.bar(x * (2 * i + 1), y, width=0.8, absolute_width=i == 1)
    return fig


@pytest.mark.mpl_image_compare
def test_bar_vectors():
    """
    Test vector arguments to bar plots.
    """
    facecolors = np.repeat(0.1, 3) * np.arange(1, 11)[:, None]
    fig, ax = pplt.subplots()
    ax.bar(
        np.arange(10),
        np.arange(1, 11),
        linewidth=3,
        edgecolor=[f"gray{i}" for i in range(9, -1, -1)],
        alpha=np.linspace(0.1, 1, 10),
        hatch=[None, "//"] * 5,
    )
    return fig


@pytest.mark.mpl_image_compare
def test_boxplot_colors():
    """
    Test box colors and cycle colors.
    """
    fig = pplt.figure(share=False)
    ax = fig.subplot(221)
    box_data = state.uniform(-3, 3, size=(1000, 5))
    violin_data = state.normal(0, 1, size=(1000, 5))
    ax.box(box_data, fillcolor=["red", "blue", "green", "orange", "yellow"])
    ax = fig.subplot(222)
    ax.violin(
        violin_data,
        fillcolor=["gray1", "gray7"],
        hatches=[None, "//", None, None, "//"],
        means=True,
        barstds=2,
    )  # noqa: E501
    ax = fig.subplot(223)
    ax.boxh(box_data, cycle="pastel2")
    ax = fig.subplot(224)
    ax.violinh(violin_data, cycle="pastel1")
    return fig


@pytest.mark.mpl_image_compare
def test_boxplot_vectors():
    """
    Test vector property arguments.
    """
    coords = (0.5, 1, 2)
    counts = (10, 20, 100)
    labels = ["foo", "bar", "baz"]
    datas = []
    for count in counts:
        data = state.rand(count)
        datas.append(data)
    datas = np.array(datas, dtype=object)
    assert len(datas) == len(coords)
    fig, ax = pplt.subplot(refwidth=3)
    ax.boxplot(
        coords,
        datas,
        lw=2,
        notch=False,
        whis=(10, 90),
        cycle="538",
        fillalpha=[0.5, 0.5, 1],
        hatch=[None, "//", "**"],
        boxlw=[2, 1, 1],
    )
    ax.format(xticklabels=labels)
    return fig


@pytest.mark.mpl_image_compare
def test_histogram_types():
    """
    Test the different histogram types using basic keywords.
    """
    fig, axs = pplt.subplots(ncols=2, nrows=2, share=False)
    data = state.normal(size=(100, 5))
    data += np.arange(5)
    kws = ({"stack": 0}, {"stack": 1}, {"fill": 0}, {"fill": 1, "alpha": 0.5})
    for ax, kw in zip(axs, kws):
        ax.hist(data, ec="k", **kw)
    return fig


@pytest.mark.mpl_image_compare
def test_invalid_plot():
    """
    Test lines with missing or invalid values.
    """
    fig, axs = pplt.subplots(ncols=2)
    data = state.normal(size=(100, 5))
    for j in range(5):
        data[:, j] = np.sort(data[:, j])
        data[: 19 * (j + 1), j] = np.nan
        # data[:20, :] = np.nan
    data_masked = ma.masked_invalid(data)  # should be same result
    for ax, dat in zip(axs, (data, data_masked)):
        ax.plot(dat, means=True, shade=True)
    return fig


@pytest.mark.mpl_image_compare
def test_invalid_dist():
    """
    Test distributions with missing or invalid data.
    """
    fig, axs = pplt.subplots(ncols=2, nrows=2)
    data = state.normal(size=(100, 5))
    for i in range(5):  # test uneven numbers of invalid values
        data[: 10 * (i + 1), :] = np.nan
    data_masked = ma.masked_invalid(data)  # should be same result
    for ax, dat in zip(axs[:2], (data, data_masked)):
        ax.violin(dat, means=True)
    for ax, dat in zip(axs[2:], (data, data_masked)):
        ax.box(dat, fill=True, means=True)
    return fig


@pytest.mark.mpl_image_compare
def test_pie_charts():
    """
    Test basic pie plots. No examples in user guide right now.
    """
    pplt.rc.inlinefmt = "svg"
    labels = ["foo", "bar", "baz", "biff", "buzz"]
    array = np.arange(1, 6)
    data = pd.Series(array, index=labels)
    fig = pplt.figure()
    ax = fig.subplot(121)
    ax.pie(array, edgefix=True, labels=labels, ec="k", cycle="reds")
    ax = fig.subplot(122)
    ax.pie(data, ec="k", cycle="blues")
    return fig


@pytest.mark.mpl_image_compare
def test_parametric_labels():
    """
    Test passing strings as parametric 'color values'. This is likely
    a common use case.
    """
    pplt.rc.inlinefmt = "svg"
    fig, ax = pplt.subplots()
    ax.parametric(
        state.rand(5), c=list("abcde"), lw=20, colorbar="b", cmap_kw={"left": 0.2}
    )
    return fig


@pytest.mark.mpl_image_compare
def test_parametric_colors():
    """
    Test color input arguments. Should be able to make monochromatic
    plots for case where we want `line` without sticky x/y edges.
    """
    fig, axs = pplt.subplots(ncols=2, nrows=2)
    colors = (
        [(0, 1, 1), (0, 1, 0), (1, 0, 0), (0, 0, 1), (1, 1, 0)],
        ["b", "r", "g", "m", "c", "y"],
        "black",
        (0.5, 0.5, 0.5),
    )
    for ax, color in zip(axs, colors):
        ax.parametric(
            state.rand(5),
            state.rand(5),
            linewidth=2,
            label="label",
            color=color,
            colorbar="b",
            legend="b",
        )
    return fig


@pytest.mark.mpl_image_compare
def test_scatter_args():
    """
    Test diverse scatter keyword parsing and RGB scaling.
    """
    x, y = state.randn(50), state.randn(50)
    data = state.rand(50, 3)
    fig, axs = pplt.subplots(ncols=4, share=0)
    ax = axs[0]
    ax.scatter(x, y, s=80, fc="none", edgecolors="r")
    ax = axs[1]
    ax.scatter(data, c=data, cmap="reds")  # column iteration
    ax = axs[2]
    with pytest.warns(pplt.internals.ProplotWarning) as record:
        ax.scatter(data[:, 0], c=data, cmap="reds")  # actual colors
    assert len(record) == 1
    ax = axs[3]
    ax.scatter(data, mean=True, shadestd=1, barstd=0.5)  # distribution
    ax.format(xlim=(-0.1, 2.1))
    return fig


@pytest.mark.mpl_image_compare
def test_scatter_inbounds():
    """
    Test in-bounds scatter plots.
    """
    fig, axs = pplt.subplots(ncols=2, share=False)
    N = 100
    fig.format(xlim=(0, 20))
    for i, ax in enumerate(axs):
        c = ax.scatter(np.arange(N), np.arange(N), c=np.arange(N), inbounds=bool(i))
        ax.colorbar(c, loc="b")
    return fig


@pytest.mark.mpl_image_compare
def test_scatter_alpha():
    """
    Test behavior with multiple alpha values.
    """
    fig, ax = pplt.subplots()
    data = state.rand(10)
    alpha = np.linspace(0.1, 1, data.size)
    ax.scatter(data, alpha=alpha)
    ax.scatter(data + 1, c=np.arange(data.size), cmap="BuRd", alpha=alpha)
    ax.scatter(data + 2, color="k", alpha=alpha)
    ax.scatter(data + 3, color=[f"red{i}" for i in range(data.size)], alpha=alpha)
    return fig


@pytest.mark.mpl_image_compare
def test_scatter_cycle():
    """
    Test scatter property cycling.
    """
    fig, ax = pplt.subplots()
    cycle = pplt.Cycle(
        "538", marker=["X", "o", "s", "d"], sizes=[20, 100], edgecolors=["r", "k"]
    )
    ax.scatter(
        state.rand(10, 4),
        state.rand(10, 4),
        cycle=cycle,
        area_size=False,
    )
    return fig


@pytest.mark.mpl_image_compare
def test_scatter_sizes():
    """
    Test marker size scaling.
    """
    # Compare setting size to input size
    size = 20
    with pplt.rc.context({"lines.markersize": size}):
        fig = pplt.figure()
        ax = fig.subplot(121, margin=0.15)
        for i in range(3):
            kw = {"absolute_size": i == 2}
            if i == 1:
                kw["smin"] = 0
                kw["smax"] = size**2  # should be same as relying on lines.markersize
            ax.scatter(np.arange(5), [0.25 * (1 + i)] * 5, size**2, **kw)
    # Test various size arguments
    ax = fig.subplot(122, margin=0.15)
    data = state.rand(5) * 500
    ax.scatter(
        np.arange(5),
        [0.25] * 5,
        c="blue7",
        sizes=[5, 10, 15, 20, 25],
        area_size=False,
        absolute_size=True,
    )
    ax.scatter(np.arange(5), [0.50] * 5, c="red7", sizes=data, absolute_size=True)
    ax.scatter(np.arange(5), [0.75] * 5, c="red7", sizes=data, absolute_size=False)
    for i, d in enumerate(data):
        ax.text(i, 0.5, format(d, ".0f"), va="center", ha="center")
    return fig
