import proplot as pplt, numpy as np, pytest


@pytest.mark.mpl_image_compare
def test_inset_basic():
    # Demonstrate that complex arrangements preserve
    # spacing, aspect ratios, and axis sharing
    gs = pplt.GridSpec(nrows=2, ncols=2)
    fig = pplt.figure(refwidth=1.5, share=False)
    for ss, side in zip(gs, "tlbr"):
        ax = fig.add_subplot(ss)
        px = ax.panel_axes(side, width="3em")
    fig.format(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="xlabel",
        ylabel="ylabel",
        xticks=0.2,
        yticks=0.2,
        title="Title",
        suptitle="Complex arrangement of panels",
        toplabels=("Column 1", "Column 2"),
        abc=True,
        abcloc="ul",
        titleloc="uc",
        titleabove=False,
    )
    return fig
