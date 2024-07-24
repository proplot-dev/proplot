#!/usr/bin/env python3
"""
Test twin, inset, and panel axes.
"""
import numpy as np

import proplot as pplt

state = np.random.RandomState(51423)


@pplt.tests.image_compare
def test_inset_colors():
    """
    Test color application for zoom boxes.
    """
    fig, ax = pplt.subplots()
    ax.format(xlim=(0, 100), ylim=(0, 100))
    ix = ax.inset_axes(
        (0.5, 0.5, 0.3, 0.3), zoom=True, zoom_kw={'color': 'r', 'fc': 'r', 'ec': 'b'}
    )  # zoom_kw={'alpha': 1})
    # ix = ax.inset_axes((40, 40, 20, 20), zoom=True, transform='data')
    ix.format(xlim=(10, 20), ylim=(10, 20), grid=False)

    fig, ax = pplt.subplots()
    ax.format(xlim=(0, 100), ylim=(0, 100))
    ix = ax.inset_axes(
        (0.3, 0.5, 0.5, 0.3), zoom=True,
        zoom_kw={'lw': 3, 'ec': 'red9', 'a': 1, 'c': pplt.set_alpha('red4', 0.5)}
    )
    ix.format(xlim=(10, 20), ylim=(10, 20))


@pplt.tests.image_compare
def test_inset_zoom_update():
    """
    Test automatic limit adjusment with successive changes. Without the extra
    lines in `draw()` and `get_tight_bbox()` this fails.
    """
    fig, ax = pplt.subplots()
    ax.format(xlim=(0, 100), ylim=(0, 100))
    ix = ax.inset_axes((40, 40, 20, 20), zoom=True, transform='data')
    ix.format(xlim=(10, 20), ylim=(10, 20), grid=False)

    ix.format(xlim=(10, 20), ylim=(10, 30))
    fig.show()

    ax.format(ylim=(0, 300))
    fig.show()


@pplt.tests.image_compare
def test_panels_with_sharing():
    """
    Previously the below text would hide the second y label.
    """
    fig, axs = pplt.subplots(ncols=2, share=False, refwidth=1.5)
    axs.panel('left')
    fig.format(ylabel='ylabel', xlabel='xlabel')


@pplt.tests.image_compare
def test_panels_without_sharing():
    """
    What should happen if `share=False` but figure-wide sharing enabled?
    Strange use case but behavior appears "correct."
    """
    fig, axs = pplt.subplots(ncols=2, share=True, refwidth=1.5, includepanels=False)
    axs.panel('left', share=False)
    fig.format(ylabel='ylabel', xlabel='xlabel')

    fig, axs = pplt.subplots(ncols=2, refwidth=1.5, includepanels=True)
    for _ in range(3):
        p = axs[0].panel('l', space=0)
        p.format(xlabel='label')
    fig.format(xlabel='xlabel')


@pplt.tests.image_compare
@pplt.tests.image_compare
def test_panels_suplabels():
    """
    Test label sharing for `includepanels=True`.
    """
    for b in range(3):
        for ncols in range(1, 3):
            fig, axs = pplt.subplots(ncols=ncols, refwidth=1.5, includepanels=(b == 2))
            if b:
                for _ in range(3):
                    axs[0].panel('l')
            axs.format(xlabel='xlabel\nxlabel\nxlabel', ylabel='ylabel', suptitle='sup')


def test_twin_axes():
    """
    Adjust twin axis positions. Should allow easily switching the location.
    """
    # Test basic twin creation and tick, spine, label location changes
    fig = pplt.figure()
    ax = fig.subplot()
    ax.format(
        ycolor='blue', ylabel='orig', ylabelcolor='blue9',
        yspineloc='l', labelweight='bold', xlabel='xlabel',
        xtickloc='t', xlabelloc='b',
    )
    ax.alty(
        loc='r', color='r', labelcolor='red9', label='other', labelweight='bold'
    )

    # Simple example but doesn't quite work. Figure out how to specify left vs. right
    # spines for 'offset' locations... maybe needs another keyword.
    fig, ax = pplt.subplots()
    ax.format(ymax=10, ylabel='Reference')
    ax.alty(color='green', label='Green', max=8)
    ax.alty(color='red', label='Red', max=15, loc=('axes', -0.2))
    ax.alty(color='blue', label='Blue', max=5, loc=('axes', 1.2), ticklabeldir='out')

    # A worked example from Riley Brady
    # Uses auto-adjusting limits
    fig, ax = pplt.subplots()
    axs = [ax, ax.twinx(), ax.twinx()]
    axs[-1].spines['right'].set_position(('axes', 1.2))
    colors = ('Green', 'Red', 'Blue')
    for ax, color in zip(axs, colors):
        data = state.random(1) * state.random(10)
        ax.plot(data, marker='o', linestyle='none', color=color)
        ax.format(ylabel='%s Thing' % color, ycolor=color)
    axs[0].format(xlabel='xlabel')
