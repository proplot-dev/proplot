#!/usr/bin/env python3
"""
Test xarray, pandas, pint, seaborn integration.
"""
import numpy as np, pandas as pd, seaborn as sns
import xarray as xr
import proplot as pplt, pytest
import pint

state = np.random.RandomState(51423)


@pytest.mark.mpl_image_compare
def test_pint_quantities():
    """
    Ensure auto-formatting and column iteration both work.
    """
    pplt.rc.unitformat = "~H"
    ureg = pint.UnitRegistry()
    fig, ax = pplt.subplots()
    ax.plot(
        np.arange(10),
        state.rand(10) * ureg.km,
        "C0",
        np.arange(10),
        state.rand(10) * ureg.m * 1e2,
        "C1",
    )
    return fig


@pytest.mark.mpl_image_compare
def test_data_keyword():
    """
    Make sure `data` keywords work properly.
    """
    N = 10
    M = 20
    ds = xr.Dataset(
        {"z": (("x", "y"), state.rand(N, M))},
        coords={
            "x": ("x", np.arange(N) * 10, {"long_name": "longitude"}),
            "y": ("y", np.arange(M) * 5, {"long_name": "latitude"}),
        },
    )
    fig, ax = pplt.subplots()
    # ax.pcolor('z', data=ds, order='F')
    ax.pcolor(z="z", data=ds, transpose=True)
    ax.format(xformatter="deglat", yformatter="deglon")
    return fig


@pytest.mark.mpl_image_compare
def test_keep_guide_labels():
    """
    Preserve metadata when passing mappables and handles to colorbar and
    legend subsequently.
    """
    fig, ax = pplt.subplots()
    df = pd.DataFrame(state.rand(5, 5))
    df.name = "variable"
    m = ax.pcolor(df)
    ax.colorbar(m)

    fig, ax = pplt.subplots()
    for k in ("foo", "bar", "baz"):
        s = pd.Series(state.rand(5), index=list("abcde"), name=k)
        ax.plot(
            s,
            legend="ul",
            legend_kw={
                "lw": 5,
                "ew": 2,
                "ec": "r",
                "fc": "w",
                "handle_kw": {"marker": "d"},
            },
        )
    return fig


@pytest.mark.mpl_image_compare
def test_seaborn_swarmplot():
    """
    Test seaborn swarm plots.
    """
    tips = sns.load_dataset("tips")
    fig = pplt.figure(refwidth=3)
    ax = fig.subplot()
    sns.swarmplot(ax=ax, x="day", y="total_bill", data=tips, palette="cubehelix")
    # fig, ax = pplt.subplots()
    # sns.swarmplot(y=state.normal(size=100), ax=ax)
    return fig


@pytest.mark.mpl_image_compare
def test_seaborn_hist():
    """
    Test seaborn histograms.
    """
    fig, axs = pplt.subplots(ncols=2, nrows=2)
    sns.histplot(state.normal(size=100), ax=axs[0])
    sns.kdeplot(x=state.rand(100), y=state.rand(100), ax=axs[1])
    penguins = sns.load_dataset("penguins")
    sns.histplot(
        data=penguins, x="flipper_length_mm", hue="species", multiple="stack", ax=axs[2]
    )
    sns.kdeplot(
        data=penguins, x="flipper_length_mm", hue="species", multiple="stack", ax=axs[3]
    )
    return fig


@pytest.mark.mpl_image_compare
def test_seaborn_relational():
    """
    Test scatter plots. Disabling seaborn detection creates mismatch between marker
    sizes and legend.
    """
    fig = pplt.figure()
    ax = fig.subplot()
    sns.set_theme(style="white")
    # Load the example mpg dataset
    mpg = sns.load_dataset("mpg")
    # Plot miles per gallon against horsepower with other semantics
    sns.scatterplot(
        x="horsepower",
        y="mpg",
        hue="origin",
        size="weight",
        sizes=(40, 400),
        alpha=0.5,
        palette="muted",
        # legend='bottom',
        # height=6,
        data=mpg,
        ax=ax,
    )
    return fig


@pytest.mark.mpl_image_compare
def test_seaborn_heatmap():
    """
    Test seaborn heatmaps. This should work thanks to backwards compatibility support.
    """
    fig, ax = pplt.subplots()
    sns.heatmap(state.normal(size=(50, 50)), ax=ax[0])
    return fig
