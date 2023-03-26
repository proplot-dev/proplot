import pytest

import numpy as np
import proplot as pplt


def test_axes_plot():
    """Test axes plots work with lists or arrays"""
    x = np.arange(10)
    fig, ax = pplt.subplots()
    ax.plot(x)
    ax.plot(x, x)

    ax.plot([1, 2, 3], [4, 5, 6])
