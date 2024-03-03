import proplot as plt, numpy as np, pytest
from matplotlib.testing import setup


@pytest.mark.mpl_image_compare
def test_standardized_inputs_1d():
    N = 5
    state = np.random.RandomState(51423)
    with plt.rc.context({"axes.prop_cycle": plt.Cycle("Grays", N=N, left=0.3)}):
        # Sample data
        x = np.linspace(-5, 5, N)
        y = state.rand(N, 5)
        fig = plt.figure(
            share=False, suptitle="Standardized input demonstration", figsize=(5, 5)
        )

        # Plot by passing both x and y coordinates
        ax = fig.subplot(121, title="Manual x coordinates")
        ax.area(x, -1 * y / N, stack=True)
        ax.bar(x, y, linewidth=0, alpha=1, width=0.8)
        ax.plot(x, y + 1, linewidth=2)
        ax.scatter(x, y + 2, marker="s", markersize=5**2)

        # Plot by passing just y coordinates
        # Default x coordinates are inferred from DataFrame,
        # inferred from DataArray, or set to np.arange(0, y.shape[0])
        ax = fig.subplot(122, title="Auto x coordinates")
        ax.area(-1 * y / N, stack=True)
        ax.bar(y, linewidth=0, alpha=1)
        ax.plot(y + 1, linewidth=2)
        ax.scatter(y + 2, marker="s", markersize=5**2)
        fig.format(xlabel="xlabel", ylabel="ylabel")
    return fig
