#!/usr/bin/env python3
"""
Test legends.
"""
import numpy as np, pandas as pd, proplot as pplt, pytest

state = np.random.RandomState(51423)


@pytest.mark.mpl_image_compare
def test_auto_legend():
    """
    Test retrieval of legends from panels, insets, etc.
    """
    fig, ax = pplt.subplots()
    ax.line(state.rand(5, 3), labels=list("abc"))
    px = ax.panel_axes("right", share=False)
    px.linex(state.rand(5, 3), labels=list("xyz"))
    # px.legend(loc='r')
    ix = ax.inset_axes((-0.2, 0.8, 0.5, 0.5), zoom=False)
    ix.line(state.rand(5, 2), labels=list("pq"))
    ax.legend(loc="b", order="F", edgecolor="red9", edgewidth=3)
    return fig


@pytest.mark.mpl_image_compare
def test_singleton_legend():
    """
    Test behavior when singleton lists are passed.
    Ensure this does not trigger centered-row legends.
    """
    fig, ax = pplt.subplots()
    h1 = ax.plot([0, 1, 2], label="a")
    h2 = ax.plot([0, 1, 1], label="b")
    ax.legend(loc="best")
    ax.legend([h1, h2], loc="bottom")
    return fig


@pytest.mark.mpl_image_compare
def test_centered_legends():
    """
    Test success of algorithm.
    """
    # Basic centered legends
    fig, axs = pplt.subplots(ncols=2, nrows=2, axwidth=2)
    hs = axs[0].plot(state.rand(10, 6))
    locs = ["l", "t", "r", "uc", "ul", "ll"]
    locs = ["l", "t", "uc", "ll"]
    labels = ["a", "bb", "ccc", "ddddd", "eeeeeeee", "fffffffff"]
    for ax, loc in zip(axs, locs):
        ax.legend(hs, loc=loc, ncol=3, labels=labels, center=True)

    # Pass centered legends with keywords or list-of-list input.
    fig, ax = pplt.subplots()
    hs = ax.plot(state.rand(10, 5), labels=list("abcde"))
    ax.legend(hs, center=True, loc="b")
    ax.legend(hs + hs[:1], loc="r", ncol=1)
    ax.legend([hs[:2], hs[2:], hs[0]], loc="t")
    return fig


@pytest.mark.mpl_image_compare
def test_manual_labels():
    """
    Test mixed auto and manual labels. Passing labels but no handles does nothing
    This is breaking change but probably best. We should not be "guessing" the
    order objects were drawn in then assigning labels to them. Similar to using
    OO interface and rejecting pyplot "current axes" and "current figure".
    """
    fig, ax = pplt.subplots()
    (h1,) = ax.plot([0, 1, 2], label="label1")
    (h2,) = ax.plot([0, 1, 1], label="label2")
    for loc in ("best", "bottom"):
        ax.legend([h1, h2], loc=loc, labels=[None, "override"])
    fig, ax = pplt.subplots()
    ax.plot([0, 1, 2])
    ax.plot([0, 1, 1])
    for loc in ("best", "bottom"):
        # ax.legend(loc=loc, labels=['a', 'b'])
        ax.legend(["a", "b"], loc=loc)  # same as above
    return fig


@pytest.mark.mpl_image_compare
def test_contour_legend_with_label():
    """
    Support contour element labels. If has no label should trigger warning.
    """
    figs = []
    label = "label"

    fig, axs = pplt.subplots(ncols=2)
    ax = axs[0]
    ax.contour(state.rand(5, 5), color="k", label=label, legend="b")
    ax = axs[1]
    ax.pcolor(state.rand(5, 5), label=label, legend="b")
    return fig


@pytest.mark.mpl_image_compare
def test_contour_legend_without_label():
    """
    Support contour element labels. If has no label should trigger warning.
    """
    figs = []
    label = None

    fig, axs = pplt.subplots(ncols=2)
    ax = axs[0]
    ax.contour(state.rand(5, 5), color="k", label=label, legend="b")
    ax = axs[1]
    ax.pcolor(state.rand(5, 5), label=label, legend="b")
    return fig


@pytest.mark.mpl_image_compare
def test_histogram_legend():
    """
    Support complex histogram legends.
    """
    pplt.rc.inlinefmt = "svg"
    fig, ax = pplt.subplots()
    res = ax.hist(
        state.rand(500, 2), 4, labels=("label", "other"), edgefix=True, legend="b"
    )
    ax.legend(res, loc="r", ncol=1)  # should issue warning after ignoring numpy arrays
    df = pd.DataFrame(
        {"length": [1.5, 0.5, 1.2, 0.9, 3], "width": [0.7, 0.2, 0.15, 0.2, 1.1]},
        index=["pig", "rabbit", "duck", "chicken", "horse"],
    )
    fig, axs = pplt.subplots(ncols=3)
    ax = axs[0]
    res = ax.hist(df, bins=3, legend=True, lw=3)
    ax.legend(loc="b")
    for ax, meth in zip(axs[1:], ("bar", "area")):
        hs = getattr(ax, meth)(df, legend="ul", lw=3)
        ax.legend(hs, loc="b")
    return fig


@pytest.mark.mpl_image_compare
def test_multiple_calls():
    """
    Test successive plotting additions to guides.
    """
    fig, ax = pplt.subplots()
    ax.pcolor(state.rand(10, 10), colorbar="b")
    ax.pcolor(state.rand(10, 5), cmap="grays", colorbar="b")
    ax.pcolor(state.rand(10, 5), cmap="grays", colorbar="b")

    fig, ax = pplt.subplots()
    data = state.rand(10, 5)
    for i in range(data.shape[1]):
        ax.plot(data[:, i], colorbar="b", label=f"x{i}", colorbar_kw={"label": "hello"})
    return fig


@pytest.mark.mpl_image_compare
def test_tuple_handles():
    """
    Test tuple legend handles.
    """
    from matplotlib import legend_handler

    fig, ax = pplt.subplots(refwidth=3, abc="A.", abcloc="ul", span=False)
    patches = ax.fill_between(state.rand(10, 3), stack=True)
    lines = ax.line(1 + 0.5 * (state.rand(10, 3) - 0.5).cumsum(axis=0))
    # ax.legend([(handles[0], lines[1])], ['joint label'], loc='bottom', queue=True)
    for hs in (lines, patches):
        ax.legend(
            [tuple(hs[:3]) if len(hs) == 3 else hs],
            ["joint label"],
            loc="bottom",
            queue=True,
            ncol=1,
            handlelength=4.5,
            handleheight=1.5,
            handler_map={tuple: legend_handler.HandlerTuple(pad=0, ndivide=3)},
        )
    return fig
