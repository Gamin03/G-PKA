"""
MODULE for NRT and MD models

@author Jimin Ma
@time   2017-6-16
"""

from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import math


def define_residual(mt, izp, iap, izt, iat):
    # compound nucleus
    kz = izt + izp
    ka = iat + iap
    izr = kz
    iar = ka
    if mt == 2:
        iar = ka - iap
        izr = kz - izp
    elif mt == 4 or (mt >= 51 and mt <= 91):
        iar = ka - 1
    elif mt == 16 or (mt >= 875 and mt <= 891):
        iar = ka - 2
    elif mt == 17:
        iar = ka - 3
    elif mt == 18:
        iar = -1
        izr = 0
    elif mt == 22:
        iar = ka - 5
        izr = kz - 2
    elif mt == 23:
        iar = ka - 13
        izr = kz - 6
    elif mt == 24:
        iar = ka - 6
        izr = kz - 2
    elif mt == 25:
        iar = ka - 7
        izr = kz - 2
    elif mt == 28:
        iar = ka - 2
        izr = kz - 1
    elif mt == 29:
        iar = ka - 9
        izr = kz - 4
    elif mt == 30:
        iar = ka - 10
        izr = kz - 4
    elif mt == 32:
        iar = ka - 3
        izr = kz - 1
    elif mt == 33:
        iar = ka - 4
        izr = kz - 1
    elif mt == 34:
        iar = ka - 4
        izr = kz - 2
    elif mt == 35:
        iar = ka - 11
        izr = kz - 5
    elif mt == 36:
        iar = ka - 12
        izr = kz - 5
    elif mt == 37:
        iar = ka - 4
    elif mt == 41:
        iar = ka - 3
        izr = kz - 1
    elif mt == 42:
        iar = ka - 4
        izr = kz - 1
    elif mt == 44:
        iar = ka - 3
        izr = kz - 2
    elif mt == 45:
        iar = ka - 6
        izr = kz - 3
    elif mt == 102:
        iar = ka
        izr = kz
    elif mt == 103 or (mt >= 600 and mt <= 649):
        iar = ka - 1
        izr = kz - 1
    elif mt == 104 or (mt >= 650 and mt <= 699):
        iar = ka - 2
        izr = kz - 1
    elif mt == 105 or (mt >= 700 and mt <= 749):
        iar = ka - 3
        izr = kz - 1
    elif mt == 106 or (mt >= 750 and mt <= 799):
        iar = ka - 3
        izr = kz - 2
    elif mt == 107 or (mt >= 800 and mt <= 849):
        iar = ka - 4
        izr = kz - 2
    elif mt == 108:
        iar = ka - 8
        izr = kz - 4
    elif mt == 109:
        iar = ka - 12
        iar = ka - 6
    elif mt == 111:
        iar = ka - 2
        izr = kz - 2
    elif mt == 112:
        iar = ka - 5
        izr = kz - 3

    return iar, izr


# Define function to calculate coefficients of NRT ane MD
#       z1, a1: projectile
#       z2, a2: material
def def_coeffs(e, mt, z1, a1, z2, a2, brk):
    constant1 = 30.724
    constant2 = 0.0793
    constant3 = 3.4008
    constant4 = 0.40244
    twothd = 2./3.
    threeq = 3./4.
    sixth = 1./6.
    onep5 = 1.5
    if e <= brk:
        return 0, 0, 0
    eL = constant1 * z1 * z2 * math.sqrt(math.pow(z1, twothd) + math.pow(z2, twothd)) * \
         (a1 + a2) / a2
    rel = 1./eL
    denom = math.pow(math.pow(z1, twothd) + math.pow(z2, twothd), threeq) * math.pow(a1, onep5) * \
            math.sqrt(a2)
    fL = constant2 * math.pow(z1, twothd) * math.sqrt(z2) * math.pow(a1 + a2, onep5) / denom
    ep = e * rel
    df = e / (1 + fL * (constant3 * math.pow(ep, sixth) + constant4 * math.pow(ep, threeq) + ep))
    eff = efficiency(z2, df, brk)
    dfeff = df * eff
    # dfeff = df * efficiency(z2, df, brk)
    return df, dfeff, float(eff)


# efficiency function for Fe from Stoller
def efficiency(z2, dam, brk):
    teff = np.array([[3.92961E+01, 3.10464E-01],
                     [4.51360E+01, 3.34810E-01],
                     [5.08440E+01, 3.47829E-01],
                     [5.64302E+01, 3.56973E-01],
                     [6.19030E+01, 3.63730E-01],
                     [1.11570E+02, 4.02940E-01],
                     [1.54243E+02, 4.22496E-01],
                     [1.91764E+02, 4.34847E-01],
                     [2.25228E+02, 4.44354E-01],
                     [2.55382E+02, 4.52342E-01],
                     [2.82773E+02, 4.55008E-01],
                     [3.07817E+02, 4.60384E-01],
                     [3.30842E+02, 4.62952E-01],
                     [3.52112E+02, 4.66695E-01],
                     [5.02267E+02, 4.76762E-01],
                     [5.91365E+02, 4.75192E-01],
                     [6.51587E+02, 4.69268E-01],
                     [6.95521E+02, 4.67951E-01],
                     [8.12090E+02, 4.62645E-01],
                     [8.97060E+02, 4.65139E-01],
                     [9.70010E+02, 4.68211E-01],
                     [1.00503E+03, 4.82998E-01],
                     [1.02928E+03, 5.01126E-01],
                     [1.05122E+03, 5.55980E-01],
                     [1.06306E+03, 6.00856E-01],
                     [1.07227E+03, 6.99076E-01],
                     [1.07676E+03, 7.77360E-01],
                     [1.07962E+03, 8.63651E-01],
                     [1.08168E+03, 9.42473E-01],
                     [1.00000E+10, 9.42473E-01]], dtype=float)
    if z2 != 26:
        return 1.0
    # for Fe-56
    emd = dam / 1000    # from eV into keV
    if emd <= 40. :
        eff = 0.5608 * math.pow(emd, -0.3029) + 3.227e-3 * emd
    else:
        f = interpolate.interp1d(teff[:, 0], teff[:, 1], kind='linear')
        eff = f(emd)

    return eff


# Find default damage displacement energy
def find_damage_displacement_energy(zz):
    edarray = {4: 31, 6: 31, 12: 25, 13: 27, 14: 25, 20: 40,
               22: 40, 23: 40, 24: 40, 25: 40, 26: 40, 27: 40, 28: 40, 29: 40,
               40: 40, 41: 40, 42: 60, 47: 60, 73: 90, 74: 90, 79: 30, 82: 25}
    if zz in edarray:
        ed = edarray[zz]
    else:
        ed = 25
    return ed