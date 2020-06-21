#!/usr/bin/env python3
"""
Tools for converting between various colorspaces. Adapted from `seaborn
<https://github.com/mwaskom/seaborn/blob/master/seaborn/external/husl.py>`__
and `hsluv-python
<https://github.com/hsluv/hsluv-python/blob/master/hsluv.py>`__.
For more info on colorspaces see the
`CIULUV specification <https://en.wikipedia.org/wiki/CIELUV>`__, the
`CIE 1931 colorspace <https://en.wikipedia.org/wiki/CIE_1931_color_space>`__,
the `HCL colorspace <https://en.wikipedia.org/wiki/HCL_color_space>`__,
and the `HSLuv system <http://www.hsluv.org/implementations/>`__.

Provided by matplotlib:

* `matplotlib.colors.to_rgb`
* `matplotlib.colors.rgb_to_hsv`
* `matplotlib.colors.hsv_to_rgb`

Provided by this module:

* `rgb_to_hsl` (same as `matplotlib.colors.rgb_to_hsv`)
* `hsl_to_rgb` (same as `matplotlib.colors.hsv_to_rgb`)
* `hcl_to_rgb`
* `rgb_to_hcl`
* `hsluv_to_rgb`
* `rgb_to_hsluv`
* `hpluv_to_rgb`
* `rgb_to_hpluv`
"""
# Imports (below functions are just meant to be used by user)
# See: https://stackoverflow.com/a/2353265/4970632
# The HLS is actually HCL
import math
from colorsys import hls_to_rgb, rgb_to_hls

# Coefficients or something
m = [
    [3.2406, -1.5372, -0.4986],
    [-0.9689, 1.8758, 0.0415],
    [0.0557, -0.2040, 1.0570]
]
m_inv = [
    [0.4124, 0.3576, 0.1805],
    [0.2126, 0.7152, 0.0722],
    [0.0193, 0.1192, 0.9505]
]
# Hard-coded D65 illuminant (has to do with expected light intensity and
# white balance that falls upon the generated color)
# See: https://en.wikipedia.org/wiki/Illuminant_D65
# Also: https://github.com/hsluv/hsluv-python/issues/3
refX = 0.95047
refY = 1.00000
refZ = 1.08883
refU = 0.19784
refV = 0.46834
lab_e = 0.008856
lab_k = 903.3


def hsluv_to_rgb(h, s, l):
    return lchuv_to_rgb(*hsluv_to_lchuv([h, s, l]))


def hsluv_to_hex(h, s, l):
    return rgb_to_hex(hsluv_to_rgb(h, s, l))


def rgb_to_hsluv(r, g, b):
    return lchuv_to_hsluv(rgb_to_lchuv(r, g, b))


def hex_to_hsluv(color):
    return rgb_to_hsluv(*hex_to_rgb(color))


def hpluv_to_rgb(h, s, l):
    return lchuv_to_rgb(*hpluv_to_lchuv([h, s, l]))


def hpluv_to_hex(h, s, l):
    return rgb_to_hex(hpluv_to_rgb(h, s, l))


def rgb_to_hpluv(r, g, b):
    return lchuv_to_hpluv(rgb_to_lchuv(r, g, b))


def hex_to_hpluv(color):
    return rgb_to_hpluv(*hex_to_rgb(color))


def lchuv_to_rgb(l, c, h):
    return CIExyz_to_rgb(CIEluv_to_CIExyz(lchuv_to_CIEluv([l, c, h])))


def rgb_to_lchuv(r, g, b):
    return CIEluv_to_lchuv(CIExyz_to_CIEluv(rgb_to_CIExyz([r, g, b])))


def hsl_to_rgb(h, s, l):
    h /= 360.0
    s /= 100.0
    l /= 100.0  # noqa
    return hls_to_rgb(h, l, s)


def rgb_to_hsl(r, g, b):
    h, l, s = rgb_to_hls(r, g, b)
    h *= 360.0
    s *= 100.0
    l *= 100.0  # noqa
    return h, s, l


def hcl_to_rgb(h, c, l):
    return CIExyz_to_rgb(CIEluv_to_CIExyz(lchuv_to_CIEluv([l, c, h])))


def rgb_to_hcl(r, g, b):
    l, c, h = CIEluv_to_lchuv(CIExyz_to_CIEluv(rgb_to_CIExyz([r, g, b])))
    return h, c, l


def rgb_prepare(triple):
    ret = []
    for ch in triple:
        ch = round(ch, 3)
        if ch < -0.0001 or ch > 1.0001:
            raise Exception(f'Illegal RGB value {ch:f}.')
        if ch < 0:
            ch = 0
        if ch > 1:
            ch = 1
        # the +0.001 fixes rounding error
        ret.append(int(round(ch * 255 + 0.001, 0)))
    return ret


def rgb_to_hex(triple):
    [r, g, b] = triple
    return '#%02x%02x%02x' % tuple(rgb_prepare([r, g, b]))


def hex_to_rgb(color):
    if color.startswith('#'):
        color = color[1:]
    r = int(color[0:2], 16) / 255.0
    g = int(color[2:4], 16) / 255.0
    b = int(color[4:6], 16) / 255.0
    return [r, g, b]


def max_chroma(L, H):
    hrad = math.radians(H)
    sinH = (math.sin(hrad))
    cosH = (math.cos(hrad))
    sub1 = (math.pow(L + 16, 3.0) / 1560896.0)
    sub2 = sub1 if sub1 > 0.008856 else (L / 903.3)
    result = float('inf')
    for row in m:
        m1 = row[0]
        m2 = row[1]
        m3 = row[2]
        top = ((0.99915 * m1 + 1.05122 * m2 + 1.14460 * m3) * sub2)
        rbottom = (0.86330 * m3 - 0.17266 * m2)
        lbottom = (0.12949 * m3 - 0.38848 * m1)
        bottom = (rbottom * sinH + lbottom * cosH) * sub2
        for t in (0.0, 1.0):
            C = (L * (top - 1.05122 * t) / (bottom + 0.17266 * sinH * t))
            if C > 0.0 and C < result:
                result = C
    return result


def hrad_extremum(L):
    lhs = (math.pow(L, 3.0) + 48.0 * math.pow(L, 2.0)
           + 768.0 * L + 4096.0) / 1560896.0
    rhs = 1107.0 / 125000.0
    sub = lhs if lhs > rhs else 10.0 * L / 9033.0
    chroma = float('inf')
    result = None
    for row in m:
        for limit in (0.0, 1.0):
            [m1, m2, m3] = row
            top = -3015466475.0 * m3 * sub + 603093295.0 * m2 * sub \
                - 603093295.0 * limit
            bottom = 1356959916.0 * m1 * sub - 452319972.0 * m3 * sub
            hrad = math.atan2(top, bottom)
            if limit == 0.0:
                hrad += math.pi
            test = max_chroma(L, math.degrees(hrad))
            if test < chroma:
                chroma = test
                result = hrad
    return result


def max_chroma_pastel(L):
    H = math.degrees(hrad_extremum(L))
    return max_chroma(L, H)


def hsluv_to_lchuv(triple):
    H, S, L = triple
    if L > 99.9999999:
        return [100, 0.0, H]
    if L < 0.00000001:
        return [0.0, 0.0, H]
    mx = max_chroma(L, H)
    C = mx * S / 100.0
    # if C > 100.0:
    #     raise ValueError(f'HSL color {triple} is outside LCH colorspace.')
    return [L, C, H]


def lchuv_to_hsluv(triple):
    L, C, H = triple
    if L > 99.9999999:
        return [H, 0.0, 100.0]
    if L < 0.00000001:
        return [H, 0.0, 0.0]
    mx = max_chroma(L, H)
    S = 100.0 * C / mx
    return [H, S, L]


def hpluv_to_lchuv(triple):
    H, S, L = triple
    if L > 99.9999999:
        return [100, 0.0, H]
    if L < 0.00000001:
        return [0.0, 0.0, H]
    mx = max_chroma_pastel(L)
    C = mx * S / 100.0
    # if C > 100.0:
    #     raise ValueError(f'HPL color {triple} is outside LCH colorspace.')
    return [L, C, H]


def lchuv_to_hpluv(triple):
    L, C, H = triple
    if L > 99.9999999:
        return [H, 0.0, 100.0]
    if L < 0.00000001:
        return [H, 0.0, 0.0]
    mx = max_chroma_pastel(L)
    S = 100.0 * C / mx
    return [H, S, L]


def dot_product(a, b):
    return sum(i * j for i, j in zip(a, b))
    # return sum(map(operator.mul, a, b))


def from_linear(c):
    if c <= 0.0031308:
        return 12.92 * c
    else:
        return (1.055 * math.pow(c, 1.0 / 2.4) - 0.055)


def to_linear(c):
    a = 0.055
    if c > 0.04045:
        return (math.pow((c + a) / (1.0 + a), 2.4))
    else:
        return (c / 12.92)


def CIExyz_to_rgb(triple):
    CIExyz = map(lambda row: dot_product(row, triple), m)
    return list(map(from_linear, CIExyz))


def rgb_to_CIExyz(triple):
    rgbl = list(map(to_linear, triple))
    return list(map(lambda row: dot_product(row, rgbl), m_inv))


def CIEluv_to_lchuv(triple):
    L, U, V = triple
    C = (math.pow(math.pow(U, 2) + math.pow(V, 2), (1.0 / 2.0)))
    hrad = (math.atan2(V, U))
    H = math.degrees(hrad)
    if H < 0.0:
        H = 360.0 + H
    return [L, C, H]


def lchuv_to_CIEluv(triple):
    L, C, H = triple
    Hrad = math.radians(H)
    U = (math.cos(Hrad) * C)
    V = (math.sin(Hrad) * C)
    return [L, U, V]


# Try setting gamma from: https://en.wikipedia.org/wiki/HCL_color_space
# The 3.0 used below should be the same; don't mess with it
gamma = 3.0  # tunable? nah, get weird stuff


def CIEfunc(t):
    if t > lab_e:
        return (math.pow(t, 1.0 / gamma))
    else:
        return (7.787 * t + 16.0 / 116.0)


def CIEfunc_inverse(t):
    if math.pow(t, 3.0) > lab_e:
        return (math.pow(t, gamma))
    else:
        return (116.0 * t - 16.0) / lab_k


def CIExyz_to_CIEluv(triple):
    X, Y, Z = triple
    if X == Y == Z == 0.0:
        return [0.0, 0.0, 0.0]
    varU = (4.0 * X) / (X + (15.0 * Y) + (3.0 * Z))
    varV = (9.0 * Y) / (X + (15.0 * Y) + (3.0 * Z))
    L = 116.0 * CIEfunc(Y / refY) - 16.0
    # Black will create a divide-by-zero error
    if L == 0.0:
        return [0.0, 0.0, 0.0]
    U = 13.0 * L * (varU - refU)
    V = 13.0 * L * (varV - refV)
    return [L, U, V]


def CIEluv_to_CIExyz(triple):
    L, U, V = triple
    if L == 0:
        return [0.0, 0.0, 0.0]
    varY = CIEfunc_inverse((L + 16.0) / 116.0)
    varU = U / (13.0 * L) + refU
    varV = V / (13.0 * L) + refV
    Y = varY * refY
    X = 0.0 - (9.0 * Y * varU) / ((varU - 4.0) * varV - varU * varV)
    Z = (9.0 * Y - (15.0 * varV * Y) - (varV * X)) / (3.0 * varV)
    return [X, Y, Z]
