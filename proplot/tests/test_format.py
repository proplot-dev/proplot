#!/usr/bin/env python3
"""
Test format and rc behavior.
"""
import locale, numpy as np, proplot as pplt, pytest
import warnings

state = np.random.RandomState(51423)


# def test_colormap_assign():
#     """
#     Test below line is possible and naming schemes.
#     """
#     pplt.rc["image.cmap"] = pplt.Colormap("phase", shift=180, left=0.2)
#     assert pplt.rc["cmap"] == pplt.rc["cmap.sequential"] == "_Phase_copy_s"
#     pplt.rc["image.cmap"] = pplt.Colormap("magma", reverse=True, right=0.8)
#     assert pplt.rc["image.cmap"] == pplt.rc["cmap.sequential"] == "_magma_copy_r"


def test_ignored_keywords():
    """
    Test ignored keywords and functions.
    """
    with warnings.catch_warnings(record=True) as record:
        fig, ax = pplt.subplots(
            gridspec_kw={"left": 3},
            subplot_kw={"proj": "cart"},
            subplotpars={"left": 0.2},
        )
    assert len(record) == 3
    with warnings.catch_warnings(record=True) as record:
        fig.subplots_adjust(left=0.2)
    assert len(record) == 1


@pytest.mark.mpl_image_compare
def test_init_format():
    """
    Test application of format args on initialization.
    """
    fig, axs = pplt.subplots(
        ncols=2,
        xlim=(0, 10),
        xlabel="xlabel",
        abc=True,
        title="Subplot title",
        collabels=["Column 1", "Column 2"],
        suptitle="Figure title",
    )
    axs[0].format(hatch="xxx", hatchcolor="k", facecolor="blue3")
    return fig


@pytest.mark.mpl_image_compare
def test_patch_format():
    """
    Test application of patch args on initialization.
    """
    fig = pplt.figure(suptitle="Super title")
    fig.subplot(
        121, proj="cyl", labels=True, land=True, latlines=20, abcloc="l", abc="[A]"
    )
    fig.subplot(
        122,
        facecolor="gray1",
        color="red",
        titleloc="l",
        title="Hello",
        abcloc="l",
        abc="[A]",
        xticks=0.1,
        yformatter="scalar",
    )
    return fig


@pytest.mark.mpl_image_compare
def test_multi_formatting():
    """
    Support formatting in multiple projections.
    """
    fig, axs = pplt.subplots(ncols=2, proj=("cart", "cyl"))
    axs[0].pcolormesh(state.rand(5, 5))
    fig.format(
        land=1,
        labels=1,
        lonlim=(0, 90),
        latlim=(0, 90),
        xlim=(0, 10),
        ylim=(0, 10),
    )
    axs[:1].format(
        land=1,
        labels=1,
        lonlim=(0, 90),
        latlim=(0, 90),
        xlim=(0, 10),
        ylim=(0, 10),
    )
    return fig


@pytest.mark.mpl_image_compare
def test_inner_title_zorder():
    """
    Test prominence of contour labels and whatnot.
    """
    fig, ax = pplt.subplots()
    ax.format(
        title="TITLE", titleloc="upper center", titleweight="bold", titlesize="xx-large"
    )
    ax.format(xlim=(0, 1), ylim=(0, 1))
    ax.text(
        0.5,
        0.95,
        "text",
        ha="center",
        va="top",
        color="red",
        weight="bold",
        size="xx-large",
    )
    x = [[0.4, 0.6]] * 2
    y = z = [[0.9, 0.9], [1.0, 1.0]]
    ax.contour(
        x,
        y,
        z,
        color="k",
        labels=True,
        levels=None,
        labels_kw={"color": "blue", "weight": "bold", "size": "xx-large"},
    )
    return fig


@pytest.mark.mpl_image_compare
def test_font_adjustments():
    """
    Test font name application. Somewhat hard to do.
    """
    fig, axs = pplt.subplots(ncols=2)
    axs.format(
        abc="A.",
        fontsize=15,
        fontname="Fira Math",
        xlabel="xlabel",
        ylabel="ylabel",
        title="Title",
        figtitle="Figure title",
        collabels=["Column 1", "Column 2"],
    )
    return fig


@pytest.mark.mpl_image_compare
def test_axes_colors():
    """
    Test behavior of passing color to format.
    """
    fig, axs = pplt.subplots(
        ncols=3,
        nrows=2,
        share=False,
        proj=("cyl", "cart", "polar", "cyl", "cart", "polar"),
        wratios=(2, 2, 1),
    )
    axs[:, 0].format(labels=True)
    axs[:3].format(edgecolor="red", gridlabelsize="med-large", gridlabelweight="bold")
    axs[:3].format(color="red")  # without this just colors the edge
    axs[1].format(xticklabelcolor="gray")
    # axs[2].format(ticklabelcolor='red')
    axs[1].format(tickcolor="blue")
    axs[3:].format(color="red")  # ensure propagates
    # axs[-1].format(gridlabelcolor='green')  # should work
    return fig


@pytest.mark.parametrize("loc", ["en_US.UTF-8"])
@pytest.mark.mpl_image_compare
def test_locale_formatting(loc):
    """
    Ensure locale formatting works. Also zerotrim should account
    for non-period decimal separators.
    """
    # dealing with read the docs
    original_locale = locale.getlocale()
    try:
        try:
            locale.setlocale(locale.LC_ALL, loc)
        except locale.Error:
            pytest.skip(f"Locale {loc} not available on this system")

        # Your test code that is sensitive to the locale settings
        assert locale.getlocale() == (loc.split(".")[0], loc.split(".")[1])

        pplt.rc["formatter.use_locale"] = False
        pplt.rc["formatter.zerotrim"] = True
        with pplt.rc.context({"formatter.use_locale": True}):
            fig, ax = pplt.subplots()
            ticks = pplt.arange(-1, 1, 0.1)
            ax.format(ylim=(min(ticks), max(ticks)), yticks=ticks)
        return fig
    finally:
        # Always reset to the original locale
        locale.setlocale(locale.LC_ALL, original_locale)
    pplt.rc["formatter.use_locale"] = False
    pplt.rc["formatter.zerotrim"] = True
    with pplt.rc.context({"formatter.use_locale": True}):
        fig, ax = pplt.subplots()
        ticks = pplt.arange(-1, 1, 0.1)
        ax.format(ylim=(min(ticks), max(ticks)), yticks=ticks)
    return fig


@pytest.mark.mpl_image_compare
def test_bounds_ticks():
    """
    Test spine bounds and location. Previously applied `fixticks`
    automatically but no longer the case.
    """
    fig, ax = pplt.subplots()
    # ax.format(xlim=(-10, 10))
    ax.format(xloc="top")
    ax.format(xlim=(-10, 15), xbounds=(0, 10))
    return fig


@pytest.mark.mpl_image_compare
def test_cutoff_ticks():
    """
    Test spine cutoff ticks.
    """
    fig, ax = pplt.subplots()
    # ax.format(xlim=(-10, 10))
    ax.format(xlim=(-10, 10), xscale=("cutoff", 0, 2), xloc="top", fixticks=True)
    ax.axvspan(0, 100, facecolor="k", alpha=0.1)
    return fig


@pytest.mark.mpl_image_compare
def test_spine_side():
    """
    Test automatic spine selection when passing `xspineloc` or `yspineloc`.
    """
    fig, ax = pplt.subplots()
    ax.plot(pplt.arange(-5, 5), (10 * state.rand(11, 5) - 5).cumsum(axis=0))
    ax.format(xloc="bottom", yloc="zero")
    ax.alty(loc="right")
    return fig


@pytest.mark.mpl_image_compare
def test_spine_offset():
    """
    Test offset axes.
    """
    fig, ax = pplt.subplots()
    ax.format(xloc="none")  # test none instead of neither
    ax.alty(loc=("axes", -0.2), color="red")
    # ax.alty(loc=('axes', 1.2), color='blue')
    ax.alty(loc=("axes", -0.4), color="blue")
    ax.alty(loc=("axes", 1.1), color="green")
    return fig


@pytest.mark.mpl_image_compare
def test_tick_direction():
    """
    Test tick direction arguments.
    """
    fig, axs = pplt.subplots(ncols=2)
    axs[0].format(tickdir="in")
    axs[1].format(xtickdirection="inout", ytickdir="out")  # rc setting should be used?
    return fig


@pytest.mark.mpl_image_compare
def test_tick_length():
    """
    Test tick length args. Ensure ratios can be applied successively.
    """
    fig, ax = pplt.subplots()
    ax.format(yticklen=100)
    ax.format(xticklen=50, yticklenratio=0.1)
    return fig


@pytest.mark.mpl_image_compare
def test_tick_width():
    """
    Test tick width args. Ensure ratios can be applied successively, setting
    width to `zero` adjusts length for label padding, and ticks can appear
    without spines if requested.
    """
    fig, axs = pplt.subplots(ncols=2, nrows=2, share=False)
    ax = axs[0]
    ax.format(linewidth=2, ticklen=20, xtickwidthratio=1)
    ax.format(ytickwidthratio=0.3)
    ax = axs[1]
    ax.format(axeslinewidth=0, ticklen=20, tickwidth=2)  # should permit ticks
    ax = axs[2]
    ax.format(tickwidth=0, ticklen=50)  # should set length to zero
    ax = axs[3]
    ax.format(linewidth=0, ticklen=20, tickwidth="5em")  # should override linewidth
    return fig


@pytest.mark.mpl_image_compare
def test_tick_labels():
    """
    Test default and overwriting properties of auto tick labels.
    """
    import pandas as pd

    data = state.rand(5, 3)
    data = pd.DataFrame(data, index=["foo", "bar", "baz", "bat", "bot"])
    fig, axs = pplt.subplots(abc="A.", abcloc="ul", ncols=2, refwidth=3, span=False)
    for i, ax in enumerate(axs):
        data.index.name = "label"
        if i == 1:
            ax.format(xformatter="null")  # overrides result
        ax.bar(data, autoformat=True)
        if i == 0:
            data.index = ["abc", "def", "ghi", "jkl", "mno"]
            data.index.name = "foobar"  # label should be updated
        ax.bar(-data, autoformat=True)
    return fig


@pytest.mark.mpl_image_compare
def test_label_settings():
    """
    Test label colors and ensure color change does not erase labels.
    """
    fig, ax = pplt.subplots()
    ax.format(xlabel="xlabel", ylabel="ylabel")
    ax.format(labelcolor="red")
    return fig
