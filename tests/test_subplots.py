import proplot as pplt, pytest


@pytest.mark.mpl_image_compare
def test_stress_test_subplots():
    fig = pplt.figure(space=0, refwidth=0.7)
    axs = fig.subplots(nrows=8, ncols=8)
    axs.format(
        abc=True,
        abcloc="ur",
        xlabel="x axis",
        ylabel="y axis",
        xticks=[],
        yticks=[],
        suptitle="A-b-c label stress test",
    )
    return fig


@pytest.mark.mpl_image_compare
def test_subplot_variable_spacing():
    # Stress test of the tight layout algorithm
    # Add large labels along the edge of one subplot
    descrip = "variable"
    equal = False
    fig, axs = pplt.subplots(
        nrows=3, ncols=3, refwidth=1.1, share=False, equal=bool(equal)
    )
    axs[1].format(xlabel="xlabel\nxlabel", ylabel="ylabel\nylabel\nylabel\nylabel")
    axs.format(
        grid=False,
        toplabels=("Column 1", "Column 2", "Column 3"),
        leftlabels=("Row 1", "Row 2", "Row 3"),
        suptitle=f"Tight layout with {descrip} row-column spacing",
    )
    return fig


@pytest.mark.mpl_image_compare
def test_subplot_equal_spacing():
    # Stress test of the tight layout algorithm
    # Add large labels along the edge of one subplot
    descrip = "equal"
    equal = True
    fig, axs = pplt.subplots(
        nrows=3, ncols=3, refwidth=1.1, share=False, equal=bool(equal)
    )
    axs[1].format(xlabel="xlabel\nxlabel", ylabel="ylabel\nylabel\nylabel\nylabel")
    axs.format(
        grid=False,
        toplabels=("Column 1", "Column 2", "Column 3"),
        leftlabels=("Row 1", "Row 2", "Row 3"),
        suptitle=f"Tight layout with {descrip} row-column spacing",
    )
    return fig
