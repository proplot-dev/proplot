#!/usr/bin/env python3
"""
Automatically build pytests from jupytext py:percent documentation files.
"""
# import glob
# import os

# import jupytext
# import pytest

# from proplot.config import rc
# from proplot.tests import SAVEFIG_KWARGS, TOLERANCE, VERSION_STRING


# def _init_tests():
#     """
#     Construct tests from the jupytext docs examples.
#     """
#     # WARNING: This will only work if all jupytext examples consist of single code
#     # cells or adjacent code cells without intervening markdown or ReST cells.
#     # Alternative would be to define the entire file as a single test but that
#     # would make testing and debugging more difficult.
#     base = os.path.dirname(__file__)
#     paths = glob.glob(f'{base}/../../docs/*.py')
#     for path in sorted(paths):
#         if os.path.basename(path) == 'conf.py':
#             continue
#         baseline, _ = os.path.splitext(os.path.basename(path))
#         result = jupytext.read(path, fmt='py:percent')
#         cells = result['cells']
#         source = ''
#         num = 1
#         for i, cell in enumerate(cells):
#             type_ = cell['cell_type']
#             if type_ == 'code':
#                 source += '\n' + cell['source']
#             if not source:
#                 continue
#             if i == len(cells) - 1 or type_ != 'code':  # end of example definition
#                 _make_cell_test(num, baseline, source)
#                 num += 1
#         print(f'\nMade {num} tests from file: {path}', end='')


# def _make_cell_test(num, baseline, cell_source):
#     """
#     Add a test using the jupytext docs cell.
#     """
#     # WARNING: Ugly kludge to replace e.g. backend='basemap' with backend='cartopy'
#     # for matplotlib versions incompatible with basemap. Generally these examples
#     # test basemap and cartopy side-by-side so this will effectively duplicate the
#     # cartopy tests but keep us from having to dump them.
#     if baseline == 'test_projections' and VERSION_STRING == 'mpl32':
#         cell_source = cell_source.replace("'basemap'", "'cartopy'")
#         if 'pplt.Proj' in cell_source or 'Table' in cell_source:
#             return  # examples that cannot be naively converted

#     def run_test():
#         rc.reset()
#         exec(cell_source)

#     name = f'{baseline}_cell_{num:02d}'
#     decorator = pytest.mark.mpl_image_compare(
#         run_test,
#         filename=f'{name}_{VERSION_STRING}.png',
#         baseline_dir='images_docs',
#         savefig_kwargs=SAVEFIG_KWARGS,
#         tolerance=TOLERANCE,
#         style={},  # no mpl style
#     )
#     test = decorator(run_test)
#     name = f'test_{name}'
#     test.__name__ = test.__qualname__ = name
#     globals()[name] = test  # for pytest detection


# # Initialize functions
# _init_tests()
