#!/usr/bin/env python3
"""
Various axis `~matplotlib.scale.ScaleBase` classes.
"""
import copy

import matplotlib.scale as mscale
import matplotlib.ticker as mticker
import matplotlib.transforms as mtransforms
import numpy as np
import numpy.ma as ma

from . import ticker as pticker
from .internals import ic  # noqa: F401
from .internals import _not_none, _version, _version_mpl, warnings

scales = mscale._scale_mapping

__all__ = [
    'CutoffScale', 'ExpScale',
    'FuncScale',
    'InverseScale',
    'LinearScale',
    'LogitScale',
    'LogScale',
    'MercatorLatitudeScale', 'PowerScale',
    'SineLatitudeScale',
    'SymmetricalLogScale',
]


def _parse_logscale_args(*keys, **kwargs):
    """
    Parse arguments for `LogScale` and `SymmetricalLogScale` that
    inexplicably require ``x`` and ``y`` suffixes by default. Also
    change the default `linthresh` to ``1``.
    """
    # NOTE: Scale classes ignore unused arguments with warnings, but matplotlib 3.3
    # version changes the keyword args. Since we can't do a try except clause, only
    # way to avoid warnings with 3.3 upgrade is to test version string.
    kwsuffix = '' if _version_mpl >= _version('3.3') else 'x'
    for key in keys:
        # Remove duplicates
        opts = {
            key: kwargs.pop(key, None),
            key + 'x': kwargs.pop(key + 'x', None),
            key + 'y': kwargs.pop(key + 'y', None),
        }
        value = _not_none(**opts)  # issues warning if multiple values passed

        # Apply defaults
        # NOTE: If linthresh is *exactly* on a power of the base, can end
        # up with additional log-locator step inside the threshold, e.g. major
        # ticks on -10, -1, -0.1, 0.1, 1, 10 for linthresh of 1. Adding slight
        # offset to *desired* linthresh prevents this.
        if key == 'linthresh' and value is None:
            value = 1 + 1e-10
        if key == 'subs' and value is None:
            value = np.arange(1, 10)
        if value is not None:  # dummy axis_name is 'x'
            kwargs[key + kwsuffix] = value

    return kwargs


class _Scale(object):
    """
    Mixin class that standardizes the behavior of
    `~matplotlib.scale.ScaleBase.set_default_locators_and_formatters`
    and `~matplotlib.scale.ScaleBase.get_transform`. Also overrides
    `__init__` so you no longer have to instantiate scales with an
    `~matplotlib.axis.Axis` instance.
    """
    def __init__(self, *args, **kwargs):
        # Pass a dummy axis to the superclass
        # WARNING: Smart bounds is deprecated. Figure out workaround by matplotlib
        # 3.4: https://github.com/matplotlib/matplotlib/pull/11004
        # Without smart bounds, inverse scale ticks disappear and Mercator ticks
        # have weird issues.
        axis = type('Axis', (object,), {'axis_name': 'x'})()
        super().__init__(axis, *args, **kwargs)
        self._default_major_locator = None
        self._default_minor_locator = None
        self._default_major_formatter = None
        self._default_minor_formatter = None

    def set_default_locators_and_formatters(self, axis, only_if_default=False):
        """
        Apply all locators and formatters defined as attributes on
        initialization and define defaults for all scales.

        Parameters
        ----------
        axis : `~matplotlib.axis.Axis`
            The axis.
        only_if_default : bool, optional
            Whether to refrain from updating the locators and formatters if
            the axis is currently using non-default versions. Useful if we
            want to avoid overwriting user customization when the scale
            is changed.
        """
        # Apply isDefault because matplotlib does this in axis._set_scale
        # but sometimes we need to bypass this method! Minor locator can be
        # "non default" even when user has not changed it, due to "turning
        # minor ticks" on and off, so set as 'default' if AutoMinorLocator.
        from .config import rc
        if not only_if_default or axis.isDefault_majloc:
            axis.set_major_locator(
                self._default_major_locator or mticker.AutoLocator()
            )
            axis.isDefault_majloc = True
        if not only_if_default or axis.isDefault_majfmt:
            axis.set_major_formatter(
                self._default_major_formatter or pticker.AutoFormatter()
            )
            axis.isDefault_majfmt = True
        if not only_if_default or axis.isDefault_minloc:
            name = axis.axis_name if axis.axis_name in 'xy' else 'x'
            axis.set_minor_locator(
                self._default_minor_locator or mticker.AutoMinorLocator()
                if rc[name + 'tick.minor.visible']
                else mticker.NullLocator()
            )
            axis.isDefault_minloc = True
        if not only_if_default or axis.isDefault_minfmt:
            axis.set_minor_formatter(
                self._default_minor_formatter or mticker.NullFormatter()
            )
            axis.isDefault_minfmt = True

    def get_transform(self):
        """
        Return the scale transform.
        """
        return self._transform


class LinearScale(_Scale, mscale.LinearScale):
    """
    As with `~matplotlib.scale.LinearScale` but with
    `~proplot.ticker.AutoFormatter` as the default major formatter.
    """
    #: The registered scale name
    name = 'linear'

    def __init__(self, **kwargs):
        """"""  # empty docstring
        super().__init__(**kwargs)
        self._transform = mtransforms.IdentityTransform()


class LogitScale(_Scale, mscale.LogitScale):
    """
    As with `~matplotlib.scale.LogitScale` but with
    `~proplot.ticker.AutoFormatter` as the default major formatter.
    """
    #: The registered scale name
    name = 'logit'

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        nonpos : {'mask', 'clip'}
          Values outside of (0, 1) can be masked as invalid, or clipped to a
          number very close to 0 or 1.
        """
        super().__init__(**kwargs)
        # self._default_major_formatter = mticker.LogitFormatter()
        self._default_major_locator = mticker.LogitLocator()
        self._default_minor_locator = mticker.LogitLocator(minor=True)


class LogScale(_Scale, mscale.LogScale):
    """
    As with `~matplotlib.scale.LogScale` but with
    `~proplot.ticker.AutoFormatter` as the default major formatter.
    ``x`` and ``y`` versions of each keyword argument are no longer
    required.
    """
    #: The registered scale name
    name = 'log'

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        base : float, optional
            The base of the logarithm. Default is ``10``.
        nonpos : {'mask', 'clip'}, optional
            Non-positive values in *x* or *y* can be masked as
            invalid, or clipped to a very small positive number.
        subs : list of int, optional
            Default *minor* tick locations are on these multiples of each power
            of the base. For example, ``subs=(1,2,5)`` draws ticks on 1, 2, 5,
            10, 20, 50, etc. The default is ``subs=numpy.arange(1, 10)``.
        basex, basey, nonposx, nonposy, subsx, subsy
            Aliases for the above keywords. These used to be conditional
            on the *name* of the axis.
        """
        keys = ('base', 'nonpos', 'subs')
        super().__init__(**_parse_logscale_args(*keys, **kwargs))
        self._default_major_locator = mticker.LogLocator(self.base)
        self._default_minor_locator = mticker.LogLocator(self.base, self.subs)


class SymmetricalLogScale(_Scale, mscale.SymmetricalLogScale):
    """
    As with `~matplotlib.scale.SymmetricalLogScale` but with
    `~proplot.ticker.AutoFormatter` as the default major formatter.
    ``x`` and ``y`` versions of each keyword argument are no longer
    required.
    """
    #: The registered scale name
    name = 'symlog'

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        base : float, optional
            The base of the logarithm. Default is ``10``.
        linthresh : float, optional
            Defines the range ``(-linthresh, linthresh)``, within which the
            plot is linear.  This avoids having the plot go to infinity around
            zero. Defaults to 2.
        linscale : float, optional
            This allows the linear range ``(-linthresh, linthresh)`` to be
            stretched relative to the logarithmic range. Its value is the
            number of decades to use for each half of the linear range. For
            example, when `linscale` is ``1`` (the default), the space used
            for the positive and negative halves of the linear range will be
            equal to one decade in the logarithmic range.
        subs : sequence of int, optional
            Default *minor* tick locations are on these multiples of each power
            of the base. For example, ``subs=(1, 2, 5)`` draws ticks on 1, 2,
            5, 10, 20, 50, 100, 200, 500, etc. The default is
            ``subs=numpy.arange(1, 10)``.
        basex, basey, linthreshx, linthreshy, linscalex, linscaley, subsx, subsy
            Aliases for the above keywords. These keywords used to be
            conditional on the name of the axis.
        """
        keys = ('base', 'linthresh', 'linscale', 'subs')
        super().__init__(**_parse_logscale_args(*keys, **kwargs))
        transform = self.get_transform()
        self._default_major_locator = mticker.SymmetricalLogLocator(transform)
        self._default_minor_locator = mticker.SymmetricalLogLocator(transform, self.subs)  # noqa: E501


class FuncScale(_Scale, mscale.ScaleBase):
    """
    An axis scale comprised of arbitrary forward and inverse transformations.
    """
    #: The registered scale name
    name = 'function'

    def __init__(
        self, arg, invert=False, parent_scale=None,
        major_locator=None, minor_locator=None,
        major_formatter=None, minor_formatter=None,
    ):
        """
        Parameters
        ----------
        arg : function, (function, function), or `~matplotlib.scale.ScaleBase`
            The transform used to translate units from the parent axis to
            the secondary axis. Input can be as follows:

            * A single function that accepts a number and returns some
              transformation of that number. If you do not provide the
              inverse, the function must be `linear \
<https://en.wikipedia.org/wiki/Linear_function>`__
              or `involutory \
<https://en.wikipedia.org/wiki/Involution_(mathematics)>`__.
              For example, to convert Kelvin to Celsius, use
              ``ax.dual%(x)s(lambda x: x - 273.15)``. To convert kilometers
              to meters, use ``ax.dual%(x)s(lambda x: x*1e3)``.
            * A 2-tuple of such functions. The second function must be the
              *inverse* of the first. For example, to apply the square, use
              ``ax.dual%(x)s((lambda x: x**2, lambda x: x**0.5))``.
              Again, if the first function is linear or involutory, you do
              not need to provide the second!
            * A scale specification interpreted by the `~proplot.constructor.Scale`
              constructor function. The forward transformation, inverse transformation,
              and default axis locators and formatters are borrowed from the resulting
              `~matplotlib.scale.ScaleBase` instance. For example, to apply the
              inverse, use ``ax.dual%(x)s('inverse')``. To apply the base-10
              exponential function, use ``ax.dual%(x)s(('exp', 10))``.

        invert : bool, optional
            If ``True``, the forward and inverse functions are *swapped*.
            Used when drawing dual axes.
        parent_scale : `~matplotlib.scale.ScaleBase`
            The axis scale of the "parent" axis. Its forward transform is
            applied to the `FuncTransform`. Used when drawing dual axes.
        major_locator, minor_locator : `~matplotlib.ticker.Locator`, optional
            The default major and minor locator. By default these are
            borrowed from `transform`. If `transform` is not an axis scale,
            they are the same as `~matplotlib.scale.LinearScale`.
        major_formatter, minor_formatter : `~matplotlib.ticker.Formatter`, optional
            The default major and minor formatter. By default these are
            borrowed from `transform`. If `transform` is not an axis scale,
            they are the same as `~matplotlib.scale.LinearScale`.
        """
        # NOTE: We permit *arbitrary* parent axis scales. If the parent is
        # non-linear, we use *its* default locators and formatters. Assumption
        # is this is a log scale and the child is maybe some multiple or offset
        # of that scale. If the parent axis scale is linear, use the funcscale
        # defaults, which can inherit defaults.
        super().__init__()
        inherit_scale = None  # scale for inheriting properties
        if callable(arg):
            forward = inverse = arg
        elif np.iterable(arg) and len(arg) == 2 and all(map(callable, arg)):
            forward, inverse = arg
        else:
            from .constructor import Scale
            try:
                inherit_scale = Scale(arg)
            except ValueError:
                raise ValueError(
                    'Input should be a function, 2-tuple of forward and and inverse '
                    f'functions, or an axis scale specification, not {arg!r}.'
                )
            trans = inherit_scale.get_transform()
            forward = trans.transform
            inverse = trans.inverted().transform

        # Create the FuncTransform or composite transform used for this class
        # May need to invert functions for dualx() and dualy()
        if invert:
            forward, inverse = inverse, forward
        functransform = FuncTransform(forward, inverse)

        # Manage the "parent" axis scale
        # NOTE: Makes sense to use the "inverse" function here because this is
        # a transformation from some *other* axis to this one, not vice versa.
        if isinstance(parent_scale, mscale.ScaleBase):
            if isinstance(parent_scale, mscale.SymmetricalLogScale):
                kwargs = {
                    key: getattr(parent_scale, key)
                    for key in ('base', 'linthresh', 'linscale', 'subs')
                }
                kwargs['linthresh'] = inverse(kwargs['linthresh'])
                parent_scale = SymmetricalLogScale(**kwargs)
            elif isinstance(parent_scale, CutoffScale):
                args = list(parent_scale.args)  # copy
                for i in range(0, len(args), 2):
                    args[i] = inverse(args[i])
                parent_scale = CutoffScale(*args)
            functransform = parent_scale.get_transform() + functransform
        elif parent_scale is not None:
            raise ValueError(
                f'parent_scale {parent_scale!r} must be a ScaleBase instance, '
                f'not {type(parent_scale)!r}.'
            )

        # Transform and default stuff
        self.functions = (forward, inverse)
        self._transform = functransform
        self._default_major_locator = major_locator
        self._default_minor_locator = minor_locator
        self._default_major_formatter = major_formatter
        self._default_minor_formatter = minor_formatter

        # Try to borrow locators and formatters
        # WARNING: Using the same locator on multiple axes can evidently
        # have unintended side effects! So we make copies.
        for scale in (inherit_scale, parent_scale):
            if not isinstance(scale, _Scale):
                continue
            if isinstance(scale, mscale.LinearScale):
                continue
            for name in (
                'major_locator', 'minor_locator', 'major_formatter', 'minor_formatter'
            ):
                attr = f'_default_{name}'
                obj = getattr(scale, attr)
                if getattr(self, attr) is None and obj is not None:
                    setattr(self, attr, copy.copy(obj))


class FuncTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, forward, inverse):
        super().__init__()
        if callable(forward) and callable(inverse):
            self._forward = forward
            self._inverse = inverse
        else:
            raise ValueError('arguments to FuncTransform must be functions')

    def inverted(self):
        return FuncTransform(self._inverse, self._forward)

    def transform_non_affine(self, values):
        with np.errstate(divide='ignore', invalid='ignore'):
            return self._forward(values)


class PowerScale(_Scale, mscale.ScaleBase):
    r"""
    "Power scale" that performs the transformation

    .. math::

        x^{c}

    """
    #: The registered scale name
    name = 'power'

    def __init__(self, power=1, inverse=False):
        """
        Parameters
        ----------
        power : float, optional
            The power :math:`c` to which :math:`x` is raised.
        inverse : bool, optional
            If ``True``, the "forward" direction performs
            the inverse operation :math:`x^{1/c}`.
        """
        super().__init__()
        if not inverse:
            self._transform = PowerTransform(power)
        else:
            self._transform = InvertedPowerTransform(power)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Return the range *vmin* and *vmax* limited to positive numbers.
        """
        if not np.isfinite(minpos):
            minpos = 1e-300
        return (
            minpos if vmin <= 0 else vmin,
            minpos if vmax <= 0 else vmax,
        )


class PowerTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, power):
        super().__init__()
        self._power = power

    def inverted(self):
        return InvertedPowerTransform(self._power)

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.power(a, self._power)


class InvertedPowerTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, power):
        super().__init__()
        self._power = power

    def inverted(self):
        return PowerTransform(self._power)

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.power(a, 1 / self._power)


class ExpScale(_Scale, mscale.ScaleBase):
    r"""
    "Exponential scale" that performs either of two transformations. When
    `inverse` is ``False`` (the default), performs the transformation

    .. math::

        Ca^{bx}

    where the constants :math:`a`, :math:`b`, and :math:`C` are set by the
    input (see below). When `inverse` is ``True``, this performs the inverse
    transformation

    .. math::

        (\log_a(x) - \log_a(C))/b

    which in appearence is equivalent to `LogScale` since it is just a linear
    transformation of the logarithm.
    """
    #: The registered scale name
    name = 'exp'

    def __init__(self, a=np.e, b=1, c=1, inverse=False):
        """
        Parameters
        ----------
        a : float, optional
            The base of the exponential, i.e. the :math:`a` in :math:`Ca^{bx}`.
        b : float, optional
            The scale for the exponent, i.e. the :math:`b` in :math:`Ca^{bx}`.
        c : float, optional
            The coefficient of the exponential, i.e. the :math:`C`
            in :math:`Ca^{bx}`.
        inverse : bool, optional
            If ``True``, the "forward" direction performs the inverse
            operation.
        """
        super().__init__()
        if not inverse:
            self._transform = ExpTransform(a, b, c)
        else:
            self._transform = InvertedExpTransform(a, b, c)

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Return the range *vmin* and *vmax* limited to positive numbers.
        """
        if not np.isfinite(minpos):
            minpos = 1e-300
        return (
            minpos if vmin <= 0 else vmin,
            minpos if vmax <= 0 else vmax,
        )


class ExpTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, a, b, c):
        super().__init__()
        self._a = a
        self._b = b
        self._c = c

    def inverted(self):
        return InvertedExpTransform(self._a, self._b, self._c)

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return self._c * np.power(self._a, self._b * np.array(a))


class InvertedExpTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, a, b, c):
        super().__init__()
        self._a = a
        self._b = b
        self._c = c

    def inverted(self):
        return ExpTransform(self._a, self._b, self._c)

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.log(a / self._c) / (self._b * np.log(self._a))


class MercatorLatitudeScale(_Scale, mscale.ScaleBase):
    """
    Axis scale that is linear in the `Mercator projection latitude \
<http://en.wikipedia.org/wiki/Mercator_projection>`__. Adapted from `this example \
<https://matplotlib.org/2.0.2/examples/api/custom_scale_example.html>`__.
    The scale function is as follows:

    .. math::

        y = \\ln(\\tan(\\pi x \\,/\\, 180) + \\sec(\\pi x \\,/\\, 180))

    The inverse scale function is as follows:

    .. math::

        x = 180\\,\\arctan(\\sinh(y)) \\,/\\, \\pi

    """
    #: The registered scale name
    name = 'mercator'

    def __init__(self, thresh=85.0):
        """
        Parameters
        ----------
        thresh : float, optional
            Threshold between 0 and 90, used to constrain axis limits between
            ``-thresh`` and ``+thresh``.
        """
        super().__init__()
        if thresh >= 90:
            raise ValueError("Mercator scale 'thresh' must be <= 90.")
        self._thresh = thresh
        self._transform = MercatorLatitudeTransform(thresh)
        self._default_major_formatter = pticker.AutoFormatter(suffix='\N{DEGREE SIGN}')  # noqa: E501

    def limit_range_for_scale(self, vmin, vmax, minpos):  # noqa: U100
        """
        Return the range *vmin* and *vmax* limited to within +/-90 degrees
        (exclusive).
        """
        return max(vmin, -self._thresh), min(vmax, self._thresh)


class MercatorLatitudeTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, thresh):
        super().__init__()
        self._thresh = thresh

    def inverted(self):
        return InvertedMercatorLatitudeTransform(self._thresh)

    def transform_non_affine(self, a):
        # NOTE: Critical to truncate valid range inside transform *and*
        # in limit_range_for_scale or get weird duplicate tick labels. This
        # is not necessary for positive-only scales because it is harder to
        # run up right against the scale boundaries.
        with np.errstate(divide='ignore', invalid='ignore'):
            m = ma.masked_where((a <= -90) | (a >= 90), a)
            if m.mask.any():
                m = np.deg2rad(m)
                return ma.log(ma.abs(ma.tan(m) + 1 / ma.cos(m)))
            else:
                a = np.deg2rad(a)
                return np.log(np.abs(np.tan(a) + 1 / np.cos(a)))


class InvertedMercatorLatitudeTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self, thresh):
        super().__init__()
        self._thresh = thresh

    def inverted(self):
        return MercatorLatitudeTransform(self._thresh)

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.rad2deg(np.arctan2(1, np.sinh(a)))


class SineLatitudeScale(_Scale, mscale.ScaleBase):
    r"""
    Axis scale that is linear in the sine transformation of *x*. The axis
    limits are constrained to fall between ``-90`` and ``+90`` degrees.
    The scale function is as follows:

    .. math::

        y = \sin(\pi x/180)

    The inverse scale function is as follows:

    .. math::

        x = 180\arcsin(y)/\pi
    """
    #: The registered scale name
    name = 'sine'

    def __init__(self):
        """"""  # no parameters
        super().__init__()
        self._transform = SineLatitudeTransform()
        self._default_major_formatter = pticker.AutoFormatter(suffix='\N{DEGREE SIGN}')

    def limit_range_for_scale(self, vmin, vmax, minpos):  # noqa: U100
        """
        Return the range *vmin* and *vmax* limited to within +/-90 degrees
        (inclusive).
        """
        return max(vmin, -90), min(vmax, 90)


class SineLatitudeTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self):
        super().__init__()

    def inverted(self):
        return InvertedSineLatitudeTransform()

    def transform_non_affine(self, a):
        # NOTE: Critical to truncate valid range inside transform *and*
        # in limit_range_for_scale or get weird duplicate tick labels. This
        # is not necessary for positive-only scales because it is harder to
        # run up right against the scale boundaries.
        with np.errstate(divide='ignore', invalid='ignore'):
            m = ma.masked_where((a < -90) | (a > 90), a)
            if m.mask.any():
                return ma.sin(np.deg2rad(m))
            else:
                return np.sin(np.deg2rad(a))


class InvertedSineLatitudeTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self):
        super().__init__()

    def inverted(self):
        return SineLatitudeTransform()

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.rad2deg(np.arcsin(a))


class CutoffScale(_Scale, mscale.ScaleBase):
    """
    Axis scale composed of arbitrary piecewise linear transformations.
    The axis can undergo discrete jumps, "accelerations", or "decelerations"
    between successive thresholds.
    """
    #: The registered scale name
    name = 'cutoff'

    def __init__(self, *args):
        """
        Parameters
        ----------
        *args : (thresh_1, scale_1, ..., thresh_N, [scale_N]), optional
            Sequence of "thresholds" and "scales". If the final scale is
            omitted (i.e. you passed an odd number of arguments) it is set
            to ``1``. Each ``scale_i`` in the sequence can be interpreted
            as follows:

            * If ``scale_i < 1``, the axis is decelerated from ``thresh_i`` to
              ``thresh_i+1``. For ``scale_N``, the axis is decelerated
              everywhere above ``thresh_N``.
            * If ``scale_i > 1``, the axis is accelerated from ``thresh_i`` to
              ``thresh_i+1``. For ``scale_N``, the axis is accelerated
              everywhere above ``thresh_N``.
            * If ``scale_i == numpy.inf``, the axis *discretely jumps* from
              ``thresh_i`` to ``thresh_i+1``. The final scale ``scale_N``
              *cannot* be ``numpy.inf``.

        Example
        -------
        >>> import proplot as pplt
        >>> import numpy as np
        >>> scale = pplt.CutoffScale(10, 0.5)  # move slower above 10
        >>> scale = pplt.CutoffScale(10, 2, 20)  # zoom out between 10 and 20
        >>> scale = pplt.CutoffScale(10, np.inf, 20)  # jump from 10 to 20
        """
        # NOTE: See https://stackoverflow.com/a/5669301/4970632
        super().__init__()
        args = list(args)
        if len(args) % 2 == 1:
            args.append(1)
        self.args = args
        self.threshs = args[::2]
        self.scales = args[1::2]
        self._transform = CutoffTransform(self.threshs, self.scales)


class CutoffTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    has_inverse = True
    is_separable = True

    def __init__(self, threshs, scales, zero_dists=None):
        # The zero_dists array is used to fill in distances where scales and
        # threshold steps are zero. Used for inverting discrete transorms.
        super().__init__()
        dists = np.diff(threshs)
        scales = np.asarray(scales)
        threshs = np.asarray(threshs)
        if len(scales) != len(threshs):
            raise ValueError(f'Got {len(threshs)} but {len(scales)} scales.')
        if any(scales < 0):
            raise ValueError('Scales must be non negative.')
        if scales[-1] in (0, np.inf):
            raise ValueError('Final scale must be finite.')
        if any(dists < 0):
            raise ValueError('Thresholds must be monotonically increasing.')
        if any((dists == 0) | (scales == 0)) and (
                any((dists == 0) != (scales == 0)) or zero_dists is None):
            raise ValueError(
                'Got zero scales and distances in different places or '
                'zero_dists is None.'
            )
        self._scales = scales
        self._threshs = threshs
        with np.errstate(divide='ignore', invalid='ignore'):
            dists = np.concatenate((threshs[:1], dists / scales[:-1]))
            if zero_dists is not None:
                dists[scales[:-1] == 0] = zero_dists
            self._dists = dists

    def inverted(self):
        # Use same algorithm for inversion!
        threshs = np.cumsum(self._dists)  # thresholds in transformed space
        with np.errstate(divide='ignore', invalid='ignore'):
            scales = 1.0 / self._scales  # new scales are inverse
        zero_dists = np.diff(self._threshs)[scales[:-1] == 0]
        return CutoffTransform(threshs, scales, zero_dists=zero_dists)

    def transform_non_affine(self, a):
        # Cannot do list comprehension because this method sometimes
        # received non-1d arrays
        dists = self._dists
        scales = self._scales
        threshs = self._threshs
        aa = np.array(a)  # copy
        with np.errstate(divide='ignore', invalid='ignore'):
            for i, ai in np.ndenumerate(a):
                j = np.searchsorted(threshs, ai)
                if j > 0:
                    aa[i] = dists[:j].sum() + (ai - threshs[j - 1]) / scales[j - 1]
        return aa


class InverseScale(_Scale, mscale.ScaleBase):
    r"""
    Axis scale that is linear in the *inverse* of *x*. The forward and inverse
    scale functions are as follows:

    .. math::

        y = x^{-1}

    """
    #: The registered scale name
    name = 'inverse'

    def __init__(self):
        """"""  # empty docstring
        super().__init__()
        self._transform = InverseTransform()
        self._default_major_locator = mticker.LogLocator(10)
        self._default_minor_locator = mticker.LogLocator(10, np.arange(1, 10))

    def limit_range_for_scale(self, vmin, vmax, minpos):
        """
        Return the range *vmin* and *vmax* limited to positive numbers.
        """
        # Unlike log-scale, we can't just warp the space between
        # the axis limits -- have to actually change axis limits. Also this
        # scale will invert and swap the limits you provide.
        if not np.isfinite(minpos):
            minpos = 1e-300
        return (
            minpos if vmin <= 0 else vmin,
            minpos if vmax <= 0 else vmax,
        )


class InverseTransform(mtransforms.Transform):
    # Create transform object
    input_dims = 1
    output_dims = 1
    is_separable = True
    has_inverse = True

    def __init__(self):
        super().__init__()

    def inverted(self):
        return InverseTransform()

    def transform_non_affine(self, a):
        with np.errstate(divide='ignore', invalid='ignore'):
            return 1.0 / a


# Monkey patch matplotlib scale factory with custom scale factory that
# accepts ScaleBase instances. This permits set_xscale and set_yscale to accept
# axis scales returned by Scale constructor and makes things constistent with
# the other constructor functions.
def _scale_factory(scale, axis, *args, **kwargs):  # noqa: U100
    """
    If `scale` is a `~matplotlib.scale.ScaleBase` instance, do nothing. If
    it is a registered scale name, look up and instantiate that scale.
    """
    if isinstance(scale, mscale.ScaleBase):
        if args or kwargs:
            warnings._warn_proplot(
                f'Ignoring args {args} and keyword args {kwargs}.'
            )
        return scale  # do nothing
    else:
        scale = scale.lower()
        if scale not in scales:
            raise ValueError(
                f'Unknown scale {scale!r}. Options are '
                + ', '.join(map(repr, scales.keys())) + '.'
            )
        return scales[scale](*args, **kwargs)


if mscale.scale_factory is not _scale_factory:
    mscale.scale_factory = _scale_factory
