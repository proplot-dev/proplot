#!/usr/bin/env python3
import numpy as np
import matplotlib.gridspec as mgridspec
default = lambda x,y: x if x is not None else y

# Generate custom GridSpec classes that override the GridSpecBase
# __setitem__ method and the 'base' __init__ method
def flexible_gridspec_factory(base):
    class _GridSpec(base):
        """
        Generalization of builtin matplotlib GridSpec that allows for
        subplots with *arbitrary spacing*. Accomplishes this by designating certain
        rows and columns as *empty*.

        Further accepts all spacing arguments in *inches*.

        This allows for user-specified extra spacing, and for automatic adjustment
        of spacing depending on whether labels or ticklabels are overlapping. Will
        be added to figure class as auto_adjust() method or something.
        """
        def __init__(self, nrows, ncols,
            hspace=None, wspace=None,   # spacing between axes
            hratios=None, wratios=None, # height and width ratios for axes
            height_ratios=None, width_ratios=None, # aliases
            **kwargs,
            ):
            # Allow passing x=None by user
            hratios = default(height_ratios, hratios)
            wratios = default(width_ratios, wratios)
            hspace = default(hspace, 0.05)
            wspace = default(wspace, 0.05)

            # Translate height/width ratios
            hratios = default(hratios, [1]*nrows)
            wratios = default(wratios, [1]*ncols)
            if len(hratios) != nrows:
                raise ValueError(f'Got {nrows} rows, but {len(hratios)} hratios.')
            if len(wratios) != ncols:
                raise ValueError(f'Got {ncols} columns, but {len(wratios)} wratios.')
            wratios = np.array(wratios)/sum(wratios)
            hratios = np.array(hratios)/sum(hratios)

            # Translate height/width spacings, implement as extra columns/rows
            try:
                if len(wspace)==1:
                    # print('hi', wspace, ncols)
                    wspace = [wspace[0]]*(ncols-1) # also convert to list
            except TypeError:
                wspace = [wspace]*(ncols-1)
            try:
                if len(hspace)==1:
                    # print('bye', hspace, nrows)
                    hspace = [hspace[0]]*(nrows-1)
            except TypeError:
                hspace = [hspace]*(nrows-1)
            if nrows>1 and len(hspace) != nrows-1:
                raise ValueError(f'Require {nrows-1} height spacings for {nrows} rows, got {len(hspace)}.')
            if ncols>1 and len(wspace) != ncols-1:
                raise ValueError(f'Require {ncols-1} width spacings for {ncols} columns, got {len(wspace)}.')
            wspace = np.array(wspace)*np.mean(wratios) # make relative to average axes width/height
            hspace = np.array(hspace)*np.mean(hratios)

            # Assign spacing as ratios
            nrows = 2*nrows - 1
            ncols = 2*ncols - 1
            wratios_final = [None]*ncols
            wratios_final[::2] = list(wratios)
            if ncols>1:
                wratios_final[1::2] = list(wspace)
            hratios_final = [None]*nrows
            hratios_final[::2] = list(hratios)
            if nrows>1:
                hratios_final[1::2] = list(hspace)

            # Continue
            return super().__init__(nrows=nrows, ncols=ncols,
                    hspace=0, wspace=0, # we implement these as invisible rows/columns
                    width_ratios=wratios_final,
                    height_ratios=hratios_final,
                    **kwargs,
                    )

        def __getitem__(self, key):
            # Interpret input
            # Note that if multiple indices are requested, e.g. gridspec[1,2], the
            # argument is passed as a tuple
            def _normalize(key, size):  # Includes last index.
                if isinstance(key, slice):
                    start, stop, _ = key.indices(size)
                    if stop > start:
                        return start, stop - 1
                else:
                    if key < 0:
                        key += size
                    if 0 <= key < size:
                        return key, key
                raise IndexError("Invalid index.")

            # Get flat (row-major) numbers
            # Note SubplotSpec initialization will figure out the row/column
            # geometry of this number automatically
            nrows, ncols = self.get_geometry()
            # print(nrows, ncols)
            if isinstance(key, tuple):
                try:
                    k1, k2 = key
                except ValueError:
                    raise ValueError('Unrecognized subplot spec "{key}".')
                num1, num2 = np.ravel_multi_index(
                    [_normalize(k1, nrows), _normalize(k2, ncols)],
                    (nrows, ncols),
                    )
                    # flags=['refs_OK'])
            else:
                num1, num2 = _normalize(key, nrows * ncols)

            # Adjust numbers to account for 'invisible' rows/columns
            def _adjust(n):
                if n<0:
                    n = nrows*ncols + n
                n = n*2 # skips the empty rows and columns
                col = n%ncols
                row = n//ncols # e.g. 3 rows and columns, index 3 is position [1,0]
                # col = col - (ncols-1)
                # row = row - (nrows-1)
                return row*ncols + col
            return mgridspec.SubplotSpec(self, _adjust(num1), _adjust(num2))
    return _GridSpec

# Make classes
FlexibleGridSpec = flexible_gridspec_factory(mgridspec.GridSpec)
FlexibleGridSpecFromSubplotSpec = flexible_gridspec_factory(mgridspec.GridSpecFromSubplotSpec)

