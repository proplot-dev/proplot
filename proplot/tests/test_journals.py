import pytest

import proplot as pplt
from proplot.subplots import JOURNAL_SPECS


# Loop through all available journals.
@pytest.mark.parametrize('journal', JOURNAL_SPECS.keys())
def test_journal_subplots(journal):
    """Tests that subplots can be generated with journal specifications."""
    f, axs = pplt.subplots(journal=journal)
