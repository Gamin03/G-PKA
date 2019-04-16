#!/usr/bin/env python
"""
Utility functions for the code.

@author Jimin Ma  <majm03@yeah.net>
@time   2018-09-11

"""
import numpy as np
import scipy.sparse as sp
from scipy import interpolate

from models import define_residual, def_coeffs, def_coeffs_njoy, find_damage_displacement_energy

# General constants
_avogadro = 6.022141930E+23
_c_light = 299792458.
_j_to_mev = 1. / 1.602176565E-13
_neutron_mass = 1.008664923
_proton_mass = 1.007825032
_np_mass_dict = {'n': 1.008664923, 'p': 1.007825032}
_np_za_dict = {'n':(0, 1), 'p': (1, 1)}


# @return (z, a) - means the daughter z and mass number tuple
def get_daughter_nuclides_particles(title, mtd, Z, A):
    particles = [('alpha', (2, 4)), ('proton', (1, 1))]
    za_incident = {'(n,': (0, 1), '(p,': (1, 1), '(a,': (2, 4)}

    # For particles
    for pt in particles:
        if pt[0] in title:
            return pt[1]

    # For recoil residual nuclide
    if "recoil" in title:
        (iz_particle, ia_particle) = za_incident[str(title[:3])]
        az_recoil = define_residual(mtd, iz_particle, ia_particle, Z, A)
    else:
        az_recoil = None

    return az_recoil[1], az_recoil[0]     # A tuple (z, a)  Attention the order!


# @return flux_pka - in pka energy structure
def interpolate_flux_pka_from_input(flux_in, flux_unit, ebound_in, ebound_pka):

    # Change the unite into per 'n s^{-1} MeV^{-1}'
    # Make sure that the interpolation happens in unit of per MeV
    flux_in_per_mev = None
    if flux_unit == 'n s^{-1}':
        flux_in_per_mev = flux_in / (ebound_in[1:] - ebound_in[:-1])
    elif flux_unit == 'n s^{-1} MeV^{-1}':
        flux_in_per_mev = flux_in

    assert (flux_in_per_mev is not None), " Flux input unit is ERROR !"

    ebound_in_mid = (ebound_in[:-1] + ebound_in[1:]) * 0.5
    # interpolate the functions
    func_interp = interpolate.interp1d(ebound_in_mid, flux_in_per_mev, kind='linear')
    ebound_pka_mid = (ebound_pka[:-1] + ebound_pka[1:]) * 0.5
    flux_pka = func_interp(ebound_pka_mid)

    # Change Unit into 'n s^{-1}', easy for matrix collapse
    flux_pka *= (ebound_pka[1:] - ebound_pka[:-1])

    return flux_pka


# @return ng_recoil_matrix - in pka energy structure
def estimate_ng_recoil_matrix(recoil_energy_group, ng_xs_array, parent_mass, key_n_p, daughter_mass):

    incident_mass = _np_mass_dict[key_n_p]

    # get the extra energy due to mass defection
    #  E = 1/2 * (\Delta m * c^2) / m_{daughter}
    #  1000. is coefficient of g to kg
    extra_energy = _j_to_mev * (parent_mass + incident_mass - daughter_mass)**2 * \
                   _c_light**2 / (2.0 * 1000.0 * _avogadro * daughter_mass)
    
    # 对每个入射能群，计算每个可能的出射能群，然后相加
    # 目前按照 SPECTER-PKA 方法来计算
    num_e_g = recoil_energy_group.shape[0] - 1
    num_ng_xs = ng_xs_array.shape[0]

    assert num_e_g == num_ng_xs, " ERROR! NG XS array size is not same as recoil energy group!"

    ng_recoil_kermas = np.zeros((num_e_g, num_e_g))
    for j in range(num_e_g):

        # 动量守恒，使用能量表述
        upper_e_out = recoil_energy_group[j+1] * incident_mass / daughter_mass
        lower_e_out = recoil_energy_group[j] * incident_mass / daughter_mass

        # 增加质量亏损引入的能量
        upper_e_out += extra_energy
        lower_e_out += extra_energy

        energy_bin_width = recoil_energy_group[1:] - recoil_energy_group[:-1]

        lower_bin = num_e_g
        upper_bin = 0
        for k in range(1, num_e_g):
            if lower_e_out < recoil_energy_group[k]:
                lower_bin = k - 1
                break

        for k in range(num_e_g - 1, 0, -1):
            if upper_e_out > recoil_energy_group[k]:
                upper_bin = k + 1
                break

        if lower_bin == num_e_g or upper_bin == 0:
            continue

        for k in range(lower_bin, upper_bin):
            ng_recoil_kermas[k, j] += ng_xs_array[j] * (min(upper_e_out, recoil_energy_group[k+1]) -
                                                        max(lower_e_out, recoil_energy_group[k]))

        # 截面xs除能量间隔
        ng_recoil_kermas[:, j] /= energy_bin_width

    M = sp.coo_matrix(ng_recoil_kermas)
    return M.T


# @parameter incident_za (0, 1) for neutron
#            target_za (Z, A) for target nuclide
# @return damage_coeffs_array - in pka energy structure
def get_damage_coeffs_array(recoil_energy_group, residual_za, target_za):
    ed = find_damage_displacement_energy(target_za[0])
    # ed = 40.0 # ====TO be removed!!!
    coeffs = []
    emid = (recoil_energy_group[1:] + recoil_energy_group[:-1]) * 0.5 * 1.E+6   # in eV
    for eg in emid[:]:
        # OLD NRT method
        # coeffs.append(def_coeffs(eg, residual_za[0], residual_za[1], target_za[0], target_za[1], ed)[0])
        # USE NJOY method
        coeffs.append(float(def_coeffs_njoy(eg, residual_za[0], residual_za[1], target_za[0], target_za[1], ed)))
    coeffs = np.array(coeffs)
    return coeffs, ed
