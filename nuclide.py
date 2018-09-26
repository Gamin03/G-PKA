#!/usr/bin/env python
"""
Class definition for nuclide.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-11

"""
import re
import string
from functools import total_ordering

import numpy as np

# Chemical element name in order of atomic number Z,
# note that element[0] is neutron 'n'
_element = 'n', \
           'H', 'He', 'Li', 'Be', 'B', \
           'C', 'N', 'O', 'F', 'Ne', \
           'Na', 'Mg', 'Al', 'Si', 'P', \
           'S', 'Cl', 'Ar', 'K', 'Ca', \
           'Sc', 'Ti', 'V', 'Cr', 'Mn', \
           'Fe', 'Co', 'Ni', 'Cu', 'Zn', \
           'Ga', 'Ge', 'As', 'Se', 'Br', \
           'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
           'Nb', 'Mo', 'Tc', 'Ru', 'Rh', \
           'Pd', 'Ag', 'Cd', 'In', 'Sn', \
           'Sb', 'Te', 'I', 'Xe', 'Cs', \
           'Ba', 'La', 'Ce', 'Pr', 'Nd', \
           'Pm', 'Sm', 'Eu', 'Gd', 'Tb', \
           'Dy', 'Ho', 'Er', 'Tm', 'Yb', \
           'Lu', 'Hf', 'Ta', 'W', 'Re', \
           'Os', 'Ir', 'Pt', 'Au', 'Hg', \
           'Tl', 'Pb', 'Bi', 'Po', 'At', \
           'Rn', 'Fr', 'Ra', 'Ac', 'Th', \
           'Pa', 'U', 'Np', 'Pu', 'Am', \
           'Cm', 'Bk', 'Cf', 'Es', 'Fm', \
           'Md', 'No', 'Lr', 'Rf', 'Db', \
           'Sg', 'Bh', 'Hs', 'Mt', 'Ds', \
           'Rg', 'Cn', 'Uut', 'Fl', 'Uup', \
           'Lv', 'Uus', 'Uuo', 'Uue', 'Ubn'

_sym2z = dict([(_element[k].upper(), k) for k in range(97)])


@total_ordering     # 让类支持比较操作
class Element:
    """
    提供元素相关信息。仅用于最终结果输出。
    """

    # Input could be 'U' or ‘92’
    def __init__(self, ele_id):
        if type(ele_id) is int:
            self.Z = ele_id
        elif type(ele_id) is str:
            ele_id = ele_id.upper()
            self.Z = _sym2z[ele_id]
        self.name = _element[self.Z]
        self.num_recoil_pka_energy_group = 0
        self.recoil_pka_energy_group = None     # np.array - Size (num_recoil_pka_energy_group + 1)
        self.recoil_pka_spectrum = None         # np.array - Size (num_recoil_pka_energy_group)
        self.average_pka_energy = 0.

        self.damage_function_coeffs = None
        self.damage_cross_section = None
        self.damage_dpa = None
        self.average_displacement_energy = 0.

    def __eq__(self, other):
        return self.Z == other.Z

    def __lt__(self, other):
        return self.Z < other.Z

    def add_recoil_pka_spectrum(self, other_recoil_pka_spectrum):
        self.recoil_pka_spectrum += other_recoil_pka_spectrum

    def copy_recoil_pka_energy_group(self, other_recoil_pka_energy_group):
        self.recoil_pka_energy_group = other_recoil_pka_energy_group

    def copy_recoil_damage(self, recoil_nuc, ratio):
        self.damage_function_coeffs = recoil_nuc.damage_function_coeffs * ratio
        self.damage_cross_section = recoil_nuc.damage_cross_section * ratio
        self.damage_dpa = recoil_nuc.damage_dpa * ratio

    def add_damage_values(self, recoil_nuc, ratio):
        self.damage_function_coeffs += recoil_nuc.damage_function_coeffs * ratio
        self.damage_cross_section += recoil_nuc.damage_cross_section * ratio
        self.damage_dpa += recoil_nuc.damage_dpa * ratio

    def calculate_average_pka_energy(self):
        if self.recoil_pka_spectrum is not None:
            data = self.recoil_pka_spectrum * (self.recoil_pka_energy_group[1:] - self.recoil_pka_energy_group[:-1])
            self.average_pka_energy = np.mean(data)


@total_ordering    # 让类支持比较操作
class Nuclide:
    """
    提供核素相关信息。仅可用于靶核、最终的核素统计两种情况。

    """

    def __init__(self, nuc_id):

        try:
            # 属性对象输入， 如 (Z, A) 或 [Z, A]
            self.Z, self.A = nuc_id.Z, nuc_id.A
        except (AttributeError, TypeError):

            try:
                # 字典类型输入，如 {'Z':92, 'A':235}
                self.Z, self.A = nuc_id['Z'], nuc_id['A']
            except (KeyError, TypeError):

                # ZAID输入： 92235
                if type(nuc_id) is int:
                    zaid = str(int(nuc_id))
                    self.Z, self.A = int(zaid[:-3]), int(zaid[-3:])

                # List输入
                if type(nuc_id) in [list, tuple]:
                    if len(nuc_id) == 2:
                        self.Z, self.A = nuc_id

                # 字符串输入
                if type(nuc_id) is str:
                    if re.search('[a-zA-Z]', nuc_id):
                        # 大写，方便比较
                        nuc_id = nuc_id.upper()

                        # 如果有连字符
                        if re.search('-', nuc_id):
                            s1, s2 = nuc_id.split('-')
                            s1 = s1.strip()
                            s2 = s2.strip()
                        else:
                            s1 = list(filter(lambda x: x in string.ascii_letters, nuc_id))
                            s2 = list(filter(lambda x: not (x in string.ascii_letters), nuc_id))
                        # 不确定s1和s2的顺序，故试一下
                        s1 = "".join(s1)
                        s2 = "".join(s2)

                        try:

                            self.Z = _sym2z[s1[:]]
                            self.A = int(s2)
                        except:
                            self.Z = _sym2z[s2[:]]
                            self.A = int(s1)

                    else:
                        zaid = str(int(nuc_id))
                        self.Z, self.A = int(zaid[:-3]), int(zaid[-3:])

        self.element = _element[self.Z]
        self.name = self.element + '-' + str(self.A)
        self.ratio = 0.
        self.mass = 0.
        self.ngamma_daughter_mass = 0.              # ONLY used in (n, gamma) matrix estimate
        self.incident_particle = 'n'                # Default particle - n
        self.num_recoil_energy_group_struc = 0
        self.recoil_energy_group_struc = None       # np.array - SIZE (num_recoil_energy_group_struc + 1)
        self.recoil_flux_pka = None                 # np.array - SIZE (num_recoil_energy_group_struc), Unit = 'n s^{-1}'

        self.recoil_nuclides_particles_info = []    # save the recoil and particle matrix info from input file
        self.recoil_nuclides = []                   # save the recoil and particle matrix
        self.ngamma_xs_array = None

        self.recoil_pka_spectrum = None             # ONLY used for total pka spectrum of this nuclide
        self.average_pka_energy = 0.

        self.damage_function_coeffs = None          # ONLY used for total damage function coeffs when needed
        self.damage_cross_section = None            # ONLY used for total damage cross sections when needed
        self.damage_dpa = None                      # ONLY used for total dpa when needed
        self.average_displacement_energy = 0.        # ONLY used for total dpa when needed

    def zaid(self):
        return self.Z * 1000 + self.A

    def __eq__(self, other):
        return (self.Z, self.A) == (other.Z, other.A)

    def __lt__(self, other):
        return (self.Z, self.A) < (other.Z, other.A)

    def set_recoil_energy_group_struc(self, eg_array):
        self.recoil_energy_group_struc = eg_array
        self.num_recoil_energy_group_struc = eg_array.shape[0] - 1

    def append_recoil_nuclide_info(self, nuc_recoil_info):
        self.recoil_nuclides_particles_info.append(nuc_recoil_info)

    def append_recoil_nuclide(self, nuc_recoil):
        self.recoil_nuclides.append(nuc_recoil)

    def set_ngamma_xs_array(self, ng_xs_array):
        self.ngamma_xs_array = ng_xs_array

    def add_recoil_pka_spectrum(self, other_recoil_pka_spectrum):
        self.recoil_pka_spectrum += other_recoil_pka_spectrum

    def copy_recoil_damage(self, recoil_nuc, ratio):
        self.damage_function_coeffs = recoil_nuc.damage_function_coeffs * ratio
        self.damage_cross_section = recoil_nuc.damage_cross_section * ratio
        self.damage_dpa = recoil_nuc.damage_dpa * ratio

    def add_damage_values(self, recoil_nuc, ratio):
        self.damage_function_coeffs += recoil_nuc.damage_function_coeffs * ratio
        self.damage_cross_section += recoil_nuc.damage_cross_section * ratio
        self.damage_dpa += recoil_nuc.damage_dpa * ratio


class NuclideRecoil:
    '''
    ONLY used for recoil nuclide of one nuclide.
    反冲核相关信息。仅用于反冲核。在最终结果统计时不能使用（转化为普通核素）。

    '''

    def __init__(self, Z, A, mtd):
        self.Z = Z
        self.A = A
        self.mtd = mtd
        self.element = _element[self.Z]
        self.name = _element[self.Z] + '-' + str(self.A)
        self.title = None
        self.mass = 0.
        self.xs_energy_group_struc = None
        self.recoil_matrix = None
        self.pka_spectrum = None

        self.estimate_ed = 0.
        self.damage_function_coeffs = None
        self.damage_cross_section = None
        self.damage_dpa = None

    # # Set the matrix energy group structure
    # def set_energy_group_structure(self, eg_array):
    #     self.xs_energy_group_struc = eg_array

    # save the recoil matrix
    def load_recoil_matrix(self, matrix):
        self.recoil_matrix = matrix

    def compute_recoil_pka_spectra(self, flux):
        B = self.recoil_matrix.toarray()
        row, column = B.shape
        row_f = flux.shape[0]

        assert (row == row_f), "PKA matrix and pka flux spectrum do not match!"

        self.pka_spectrum = np.dot(B.T, flux)


if __name__ == '__main__':

    class Foo:
       pass

    nuc_obj = Foo()
    nuc_obj.Z = 92
    nuc_obj.A = 235

    nuc_ids =   [ 'U235', 'U-235', '235U', '235-U',
                  'u235', 'u-235', '235u', '235-u',
                   92235, "92235",
                   (92,235), [92, 235],
                   {'Z':92, 'A':235},
                   nuc_obj
                ]

    for nuc_id in nuc_ids:
        nuclide = Nuclide(nuc_id)
        print(nuc_id, type(nuc_id), nuclide.Z, nuclide.A, nuclide.element, nuclide.name)



        assert nuclide.Z == 92
        assert nuclide.A == 235
        assert nuclide.element == 'U'