#!/usr/bin/env python3
"""
Test twin, inset, and panel axes.
"""
import numpy as np
import pytest
import proplot as pplt

state = np.random.RandomState(51423)


def test_axis_access():
    # attempt to access the ax object 2d and linearly
    fix, ax = pplt.subplots(ncols=2, nrows=2)
    ax[0, 0]
    ax[1, 0]
    with pytest.raises(IndexError):
        ax[0, 3]
    ax[3]


@pytest.mark.mpl_image_compare
def test_inset_colors_1():
    """
    Test color application for zoom boxes.
    """
    fig, ax = pplt.subplots()
    ax.format(xlim=(0, 100), ylim=(0, 100))
    ix = ax.inset_axes(
        (0.5, 0.5, 0.3, 0.3), zoom=True, zoom_kw={"color": "r", "fc": "r", "ec": "b"}
    )  # zoom_kw={'alpha': 1})
    # ix = ax.inset_axes((40, 40, 20, 20), zoom=True, transform='data')
    ix.format(xlim=(10, 20), ylim=(10, 20), grid=False)
    return fig


@pytest.mark.mpl_image_compare
def test_inset_colors_2():
    fig, ax = pplt.subplots()
    ax.format(xlim=(0, 100), ylim=(0, 100))
    ix = ax.inset_axes(
        (0.3, 0.5, 0.5, 0.3),
        zoom=True,
        zoom_kw={"lw": 3, "ec": "red9", "a": 1, "c": pplt.set_alpha("red4", 0.5)},
    )
    ix.format(xlim=(10, 20), ylim=(10, 20))
    return fig


@pytest.mark.mpl_image_compare
def test_inset_zoom_update():
    """
    Test automatic limit adjusment with successive changes. Without the extra
    lines in `draw()` and `get_tight_bbox()` this fails.
    """
    fig, ax = pplt.subplots()
    ax.format(xlim=(0, 100), ylim=(0, 100))
    ix = ax.inset_axes((40, 40, 20, 20), zoom=True, transform="data")
    ix.format(xlim=(10, 20), ylim=(10, 20), grid=False)
    ix.format(xlim=(10, 20), ylim=(10, 30))
    ax.format(ylim=(0, 300))
    return fig


@pytest.mark.mpl_image_compare
def test_panels_with_sharing():
    """
    Previously the below text would hide the second y label.
    """
    fig, axs = pplt.subplots(ncols=2, share=False, refwidth=1.5)
    axs.panel("left")
    fig.format(ylabel="ylabel", xlabel="xlabel")
    return fig


@pytest.mark.mpl_image_compare
def test_panels_without_sharing_1():
    """
    What should happen if `share=False` but figure-wide sharing enabled?
    Strange use case but behavior appears "correct."
    """
    fig, axs = pplt.subplots(ncols=2, share=True, refwidth=1.5, includepanels=False)
    axs.panel("left", share=False)
    fig.format(ylabel="ylabel", xlabel="xlabel")
    return fig


@pytest.mark.mpl_image_compare
def test_panels_without_sharing_2():
    fig, axs = pplt.subplots(ncols=2, refwidth=1.5, includepanels=True)
    for _ in range(3):
        p = axs[0].panel("l", space=0)
        p.format(xlabel="label")
    fig.format(xlabel="xlabel")
    return fig


@pytest.mark.mpl_image_compare
def test_panels_suplabels_three_hor_panels():
    """
    Test label sharing for `includepanels=True`.
    Test for 1 subplot with 3 left panels
    Include here centers the x label to include the panels
    The xlabel should be centered along the main plot with the included side panels
    """
    fig = pplt.figure()
    ax = fig.subplots(refwidth=1.5, includepanels=True)
    for _ in range(3):
        ax[0].panel("l")
    ax.format(xlabel="xlabel", ylabel="ylabel\nylabel\nylabel", suptitle="sup")
    return fig


def test_panels_suplabels_three_hor_panels_donotinlcude():
    """
    Test label sharing for `includepanels=True`.
    Test for 1 subplot with 3 left panels
    The xlabel should be centered on the main plot
    """
    fig = pplt.figure()
    ax = fig.subplots(refwidth=1.5, includepanels=False)
    for _ in range(3):
        ax[0].panel("l")
    ax.format(
        xlabel="xlabel",
        ylabel="ylabel\nylabel\nylabel",
        suptitle="sup",
    )
    return fig


@pytest.mark.mpl_image_compare
def test_twin_axes_1():
    """
    Adjust twin axis positions. Should allow easily switching the location.
    """
    # Test basic twin creation and tick, spine, label location changes
    fig = pplt.figure()
    ax = fig.subplot()
    ax.format(
        ycolor="blue",
        ylabel="orig",
        ylabelcolor="blue9",
        yspineloc="l",
        labelweight="bold",
        xlabel="xlabel",
        xtickloc="t",
        xlabelloc="b",
    )
    ax.alty(loc="r", color="r", labelcolor="red9", label="other", labelweight="bold")
    return fig


@pytest.mark.mpl_image_compare
def test_twin_axes_2():
    # Simple example but doesn't quite work. Figure out how to specify left vs. right
    # spines for 'offset' locations... maybe needs another keyword.
    fig, ax = pplt.subplots()
    ax.format(ymax=10, ylabel="Reference")
    ax.alty(color="green", label="Green", max=8)
    ax.alty(color="red", label="Red", max=15, loc=("axes", -0.2))
    ax.alty(color="blue", label="Blue", max=5, loc=("axes", 1.2), ticklabeldir="out")
    return fig


@pytest.mark.mpl_image_compare
def test_twin_axes_3():
    # A worked example from Riley Brady
    # Uses auto-adjusting limits
    fig, ax = pplt.subplots()
    axs = [ax, ax.twinx(), ax.twinx()]
    axs[-1].spines["right"].set_position(("axes", 1.2))
    colors = ("Green", "Red", "Blue")
    for ax, color in zip(axs, colors):
        data = state.random(1) * state.random(10)
        ax.plot(data, marker="o", linestyle="none", color=color)
        ax.format(ylabel="%s Thing" % color, ycolor=color)
    axs[0].format(xlabel="xlabel")
    return fig
