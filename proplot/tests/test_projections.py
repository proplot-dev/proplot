#!/usr/bin/env python3
"""
Test projection features.
"""
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import proplot as pplt
import pytest

state = np.random.RandomState(51423)


@pytest.mark.mpl_image_compare
def test_aspect_ratios():
    """
    Test aspect ratio adjustments.
    """
    fig, axs = pplt.subplots(ncols=2)
    axs[0].format(aspect=1.5)
    fig, axs = pplt.subplots(ncols=2, proj=("cart", "cyl"), aspect=2)
    axs[0].set_aspect(1)
    return fig


if pplt.internals._version_mpl <= "3.2":

    @pytest.mark.mpl_image_compare
    def test_basemap_labels():
        """
        Add basemap labels.
        """
        fig, axs = pplt.subplots(ncols=2, proj="robin", refwidth=3, basemap=True)
        axs.format(coast=True, labels="rt")
        return fig


@pytest.mark.mpl_image_compare
def test_cartopy_labels():
    """
    Add cartopy labels.
    """
    fig, axs = pplt.subplots(ncols=2, proj="robin", refwidth=3)
    axs.format(coast=True, labels=True)
    axs[0].format(inlinelabels=True)
    axs[1].format(rotatelabels=True)
    return fig


@pytest.mark.mpl_image_compare
def test_cartopy_contours():
    """
    Test bug with cartopy contours. Sometimes results in global coverage
    with single color sometimes not.
    """
    N = 10
    fig = plt.figure(figsize=(5, 2.5))
    ax = fig.add_subplot(projection=ccrs.Mollweide())
    ax.coastlines()
    x = np.linspace(-180, 180, N)
    y = np.linspace(-90, 90, N)
    z = state.rand(N, N) * 10 - 5
    m = ax.contourf(
        x,
        y,
        z,
        transform=ccrs.PlateCarree(),
        cmap="RdBu_r",
        vmin=-5,
        vmax=5,
    )
    fig.colorbar(m, ax=ax)
    fig = pplt.figure()
    ax = fig.add_subplot(projection=pplt.Mollweide(), extent="auto")
    ax.coastlines()
    N = 10
    m = ax.contourf(
        np.linspace(0, 180, N),
        np.linspace(0, 90, N)[1::2],
        state.rand(N // 2, N) * 10 + 5,
        cmap="BuRd",
        transform=pplt.PlateCarree(),
        edgefix=False,
    )
    fig.colorbar(m, ax=ax)
    return fig


@pytest.mark.mpl_image_compare
def test_cartopy_manual():
    """
    Test alternative workflow without classes.
    """
    fig = pplt.figure()
    proj = pplt.Proj("npstere")
    # fig.add_axes([0.1, 0.1, 0.9, 0.9], proj='geo', map_projection=proj)
    fig.add_subplot(111, proj="geo", land=True, map_projection=proj)
    return fig


@pytest.mark.mpl_image_compare
def test_three_axes():
    """
    Test basic 3D axes here.
    """
    pplt.rc["tick.minor"] = False
    fig, ax = pplt.subplots(proj="3d", outerpad=3)
    return fig


@pytest.mark.mpl_image_compare
def test_projection_dicts():
    """
    Test projection dictionaries.
    """
    fig = pplt.figure(refnum=1)
    a = [[1, 0], [1, 4], [2, 4], [2, 4], [3, 4], [3, 0]]
    fig.subplots(a, proj={1: "cyl", 2: "cart", 3: "cart", 4: "cart"})
    return fig


@pytest.mark.mpl_image_compare
def test_polar_projections():
    """
    Rigorously test polar features here.
    """
    fig, ax = pplt.subplots(proj="polar")
    ax.format(
        rlabelpos=45,
        thetadir=-1,
        thetalines=90,
        thetalim=(0, 270),
        theta0="N",
        r0=0,
        rlim=(0.5, 1),
        rlines=0.25,
    )
    return fig
