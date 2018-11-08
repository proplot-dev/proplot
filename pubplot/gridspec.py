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
        def __init__(self, nrows, ncols, **kwargs):
            # Add these as attributes; want _ratios to be self-contained, so
            # it can be invoked on already instantiated gridspec
            self._nrows_visible = nrows
            self._ncols_visible = ncols
            self._nrows = nrows*2-1
            self._ncols = ncols*2-1
            # print('before', kwargs)
            wratios, hratios, kwargs = self._ratios(**kwargs)
            # print(self._nrows, self._ncols, len(wratios), len(hratios))
            # print('after', kwargs)
            return super().__init__(self._nrows, self._ncols,
                    hspace=0, wspace=0, # we implement these as invisible rows/columns
                    width_ratios=wratios,
                    height_ratios=hratios,
                    **kwargs,
                    )

        def __getitem__(self, key):
            # Interpret input
            # Note that if multiple indices are requested, e.g. gridspec[1,2], the
            # argument is passed as a tuple
            def _normalize(key, size):
                if isinstance(key, slice):
                    start, stop, _ = key.indices(size)
                    if stop > start:
                        return start, stop - 1
                else:
                    if key < 0:
                        key += size
                    if 0 <= key < size:
                        return key, key
                raise IndexError(f"Invalid index: {key} with size {size}.")
            # Get flat (row-major) numbers
            # Note SubplotSpec initialization will figure out the row/column
            # geometry of this number automatically
            # nrows, ncols = self.get_geometry()
            # nrows, ncols = self._nrows_visible, self._ncols_visible
            nrows, ncols = self._nrows, self._ncols
            nrows_visible, ncols_visible = self._nrows_visible, self._ncols_visible
            if isinstance(key, tuple):
                try:
                    k1, k2 = key
                except ValueError:
                    raise ValueError('Unrecognized subplot spec "{key}".')
                num1, num2 = np.ravel_multi_index(
                    [_normalize(k1, nrows_visible), _normalize(k2, ncols_visible)],
                    (nrows, ncols),
                    )
            else:
                num1, num2 = _normalize(key, nrows_visible * ncols_visible)
            # Adjust numbers to account for 'invisible' rows/columns
            # Think proceeding down by rows, when you move to a new column
            # that skips a 'hspace' and when you move to a new row that skips
            # a 'wspace', so just multiply by 2!
            def _adjust(n):
                if n<0:
                    n = 2*(n+1) - 1 # want -1 to stay -1, -2 becomes -3, etc.
                else:
                    n = n*2
                return n
            # print('nums', nrows, ncols, num1, num2, _adjust(num1), _adjust(num2))
            num1, num2 = _adjust(num1), _adjust(num2)
            return mgridspec.SubplotSpec(self, num1, num2)

        def _ratios(self, hspace=None, wspace=None, # spacing between axes
                hratios=None, wratios=None,
                height_ratios=None, width_ratios=None,
                **kwargs):
            # Allow passing x=None by user
            nrows = self._nrows_visible
            ncols = self._ncols_visible
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
            wratios_final = [None]*self._ncols
            wratios_final[::2] = list(wratios)
            if self._ncols>1:
                wratios_final[1::2] = list(wspace)
            hratios_final = [None]*self._nrows
            hratios_final[::2] = list(hratios)
            if self._nrows>1:
                hratios_final[1::2] = list(hspace)
            return wratios_final, hratios_final, kwargs # bring extra kwargs back

        def update(self, **kwargs):
            # Handle special hspace/wspace arguments, and just set the simple
            # left/right/top/bottom attributes
            wratios, hratios, kwargs = self._ratios(**kwargs)
            self.set_width_ratios(wratios)
            self.set_height_ratios(hratios)
            for key,value in kwargs.items():
                setattr(self,key,value)
            return super().update()

    return _GridSpec

# Make classes
FlexibleGridSpec = flexible_gridspec_factory(mgridspec.GridSpec)
FlexibleGridSpecFromSubplotSpec = flexible_gridspec_factory(mgridspec.GridSpecFromSubplotSpec)

