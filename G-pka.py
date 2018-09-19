#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main for G-pka.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-11
@note   Rewritten the old G-PKA.
        - Add class definitions
        - Add the input file.
        - Change the pka xs matrix file into SPECTER-PKA format.

@warning
        - [未解决的问题]
        2018/09/14 1) 对于算例，元素Zr的结果，PKA谱型相似，但比SPECTER-PKA要大1-2个量级，有待查找原因！
                      => 已解决。原因是注意输出是否皈依到一个核素。
                   2) PKA能谱需要归一到一个靶原子，然后进行输出
                      => 已解决。2018/09/16
                   3) (n,gamma)待处理
                      => 已解决。采用SPECTER-PKA完全相同的方法。
                   4) Zr-93,95,97 共振区域结果不正确！
                      => 已解决。注意(n,gamma)散射矩阵的入射和出射能群下标，以及计算特定群的上下群边界值，这些容易出错。
        2018/09/19 1) 增加dpa损伤计算功能
                   2) 输出结果格式的选择，输出最终PKA和dpa结果数据
                   3) 输出绘图格式，绘图的多种形式，绘图到文件

"""
# import sys
# import numpy as np
# import matplotlib.pyplot as plt

from nuclide import Element, Nuclide, NuclideRecoil
from input import Input
from read_pka_file import *
from utility_pka import *
from utility_fig import *
from utility_basic import *

# Code version control
_major_version = 0
_minor_version = 1
_bugfix_version = 0
print_code_title_version(_major_version, _minor_version, _bugfix_version)

# Global variables
global_recoil = {}
global_element = {}



# --- Initialization ----------------------------------------------------------

if len(sys.argv) == 1:
    input_file = 'input.json'
    output_file = 'output.txt'
elif len(sys.argv) == 2:
    input_file = str(sys.argv[0])
    output_file = 'output.txt'
elif len(sys.argv) == 3:
    input_file = str(sys.argv[0])
    output_file = str(sys.argv[1])
else:
    raise ValueError('Input parameters ERROR!')

inp = Input(input_file)
output = open(output_file, 'w')
nuclides = []

# Print title into output
print_code_title_version(_major_version, _minor_version, _bugfix_version, file=output)

energy_group = None  # =====Temp

# --- Read the input file -----------------------------------------------------

inp.read_infile(output)
inp.read_flux(output)

# --- Read the pka xs matrix files --------------------------------------------
print("\n\n>>> START READ PKA MATRIX FILES ...", file=output)
for i in range(inp.number_pka_files):
    # Set up a new parent nuclide and its ratio
    nuc = Nuclide(inp.pka_files[i]['parent'])
    nuc.ratio = inp.pka_files[i]['pka_ratios']

    # Read the parent and daughter mass for ngamma reaction channel
    nuc.mass = inp.pka_files[i]['ngamma_parent_mass']
    nuc.ngamma_daughter_mass = inp.pka_files[i]['ngamma_daughter_mass']

    print("\tPKA FILE    : {}".format(inp.pka_files[i]['pka_filename']), file=output)
    print("\tNUCLIDE INFO: {} Z={}, A={}, ratio={}".format(nuc.name, nuc.Z, nuc.A, nuc.ratio), file=output)

    # Read the parent nuclide pka file
    with open(inp.pka_files[i]['pka_filename']) as pka_file:

        print("    >>> LIST OF MTD WHICH HAS BEEN READ:", file=output)

        # 读取第一部分，能量网格结构
        title, mtd, pi, pe_array, A = read_pka_file_energy_group_struc(pka_file)
        nuc.num_recoil_energy_group_struc = pi
        nuc.set_recoil_energy_group_struc(pe_array)
        if energy_group is None:
            energy_group = pe_array
        nuc.append_recoil_nuclide_info((title, mtd, A))   # Add the recoil info
        print("\t\t| {0:3d} | {1:30s} |".format(mtd, title), file=output)

        # 读取第二部分：多个不同的反冲核、生成粒子的群群矩阵，稀疏矩阵格式
        while True:
            title, mtd, A = read_pka_file_each_matrix(pka_file)
            if A is None:
                break
            else:
                nuc.append_recoil_nuclide_info((title, mtd, A))
                print("\t\t| {0:3d} | {1:30s} |".format(mtd, title), file=output)

        # 读取第三部分：(n, g) 反应截面。单独存储。
        title, mtd, ng_xs_array = read_pka_file_ng_xs(pka_file)

        # 若需要处理 (n,g) 反应，则在此处生成 (n,g) 群群矩阵
        if inp.do_gamma_estimate and mtd == 102 and 'cross' in title:
            nuc.set_ngamma_xs_array(ng_xs_array)
            A = estimate_ng_recoil_matrix(nuc.recoil_energy_group_struc, nuc.ngamma_xs_array, nuc.mass,
                                          nuc.incident_particle, nuc.ngamma_daughter_mass)
            nuc.append_recoil_nuclide_info(('(n,g) recoil matrix', mtd, A))
            print("\t\t| {0:3d} | {1:30s} |".format(mtd, '(n,g) recoil matrix [estimated]'), file=output)

    nuclides.append(nuc)

# --- Calculate PKA and DPA values --------------------------------------------
print("\n\n>>> START CALCULATE PKA AND DPA VALUES ...", file=output)
for nuc in nuclides:
    # prepare for collapse
    nuc.recoil_flux_pka = interpolate_flux_pka_from_input(inp.flux_spectrum,
                                inp.flux_unit,
                                inp.flux_energy_group,
                                nuc.recoil_energy_group_struc)
    recoil_pka_spectrum_total = np.zeros(nuc.num_recoil_energy_group_struc)

    # Deal each nuclide
    for k, recoil_nuc in enumerate(nuc.recoil_nuclides_particles_info):

        # Deal each reaction channel
        # --- Set up the NuclideRecoil for this nuc
        za_recoil = get_daughter_nuclides_particles(recoil_nuc[0], recoil_nuc[1], nuc.Z, nuc.A)
        nuc_recoil = NuclideRecoil(za_recoil[0], za_recoil[1], recoil_nuc[1])
        # --- Load the recoil matrix
        nuc_recoil.load_recoil_matrix(recoil_nuc[2])
        # --- Calculate recoil pka spectra and save it
        nuc_recoil.compute_recoil_pka_spectra(nuc.recoil_flux_pka)
        nuc.append_recoil_nuclide(nuc_recoil)

        # Accumulate the NuclideRecoil for this nuc
        recoil_pka_spectrum_total += nuc_recoil.pka_spectrum

    # Save the total recoil pka spectrum for this nuc
    nuc.recoil_pka_spectrum = recoil_pka_spectrum_total

# --- Output results ----------------------------------------------------------
# Check each recoil nuclides of every initial nuclide, save them into global
# nuclides and elements dictionary
# ****** Add the nuclide ratio here ******
print("\n\n>>> START CALCULATE TOTAL PKA RESULTS ...", file=output)
for nuc in nuclides:
    ratio = nuc.ratio
    for recoil_nuc in nuc.recoil_nuclides:

        # Check if the nuclide exist in the dict, if not, append it
        if recoil_nuc.name in global_recoil:
            if ( (recoil_nuc.name != 'He-4') and (800 <= recoil_nuc.mtd <= 849) ) or \
               ( (recoil_nuc.name != 'H-1')  and (600 <= recoil_nuc.mtd <= 649) ) :
               continue
            global_recoil[recoil_nuc.name].add_recoil_pka_spectrum(recoil_nuc.pka_spectrum * ratio)
        else:
            # Increase new nuclide into global nuclide dictionary
            nuclide = Nuclide((recoil_nuc.Z, recoil_nuc.A))
            nuclide.recoil_pka_spectrum = recoil_nuc.pka_spectrum * ratio
            global_recoil[recoil_nuc.name] = nuclide

        # Check if the element exist in the dict, if not, append it
        if recoil_nuc.element in global_element:
            if ( (recoil_nuc.name != 'He-4') and (800 <= recoil_nuc.mtd <= 849) ) or \
               ( (recoil_nuc.name != 'H-1')  and (600 <= recoil_nuc.mtd <= 649) ) :
               continue
            global_element[recoil_nuc.element].add_recoil_pka_spectrum(recoil_nuc.pka_spectrum * ratio)
        else:
            # Increase new element into global element dictionary
            element = Element(recoil_nuc.Z)
            element.num_recoil_pka_energy_group = nuc.num_recoil_energy_group_struc
            element.copy_recoil_pka_energy_group(nuc.recoil_energy_group_struc)
            element.recoil_pka_spectrum = recoil_nuc.pka_spectrum * ratio
            global_element[element.name] = element

# ------ Sum results ------
total_pka_spectrum = None
for key in global_element:
    ele = global_element[key]
    ele.calculate_average_pka_energy()
    if total_pka_spectrum is None:
        total_pka_spectrum = ele.recoil_pka_spectrum
    else:
        total_pka_spectrum += ele.recoil_pka_spectrum


print("\n\n>>> OUTPUT TOTAL PKA RESULTS ...", file=output)

element = global_element['Zr']
print("\n#### {} ".format(element.name), file=output)
for i in range(element.num_recoil_pka_energy_group):
    print(" {0:12.4E} {1:12.4E} {2:12.4E}".format(element.recoil_pka_energy_group[i],
                                                  element.recoil_pka_energy_group[i+1],
                                                  element.recoil_pka_spectrum[i]), file=output)
print("#    Average PKA energy = {}".format(element.average_pka_energy), file=output)
output.close()


density = 6.506
atomic_mass = 91.22364157600980
atoms_per_mole = 0.6022

# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# # ax.loglog(energy_group[:-1], global_element['Zr'].recoil_pka_spectrum)
# ax.loglog(energy_group[:-1], global_recoil['H-1'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='g', label='H-1')
# ax.loglog(energy_group[:-1], global_recoil['He-4'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='m', label='He-4')
# ax.loglog(energy_group[:-1], global_recoil['Zr-90'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='b', label='Zr-90')
# ax.loglog(energy_group[:-1], global_recoil['Zr-91'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='k', label='Zr-91')
# ax.loglog(energy_group[:-1], global_recoil['Zr-92'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='y', label='Zr-92')
# ax.loglog(energy_group[:-1], global_recoil['Zr-94'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='c', label='Zr-94')
# ax.loglog(energy_group[:-1], global_recoil['Zr-95'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='b', label='Zr-95')
# ax.loglog(energy_group[:-1], global_recoil['Zr-96'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='m', label='Zr-96')
# ax.loglog(energy_group[:-1], global_recoil['Zr-97'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='b', label='Zr-97')
# ax.loglog(energy_group[:-1], global_recoil['Zr-93'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='g', label='Zr-93')
# ax.loglog(energy_group[:-1], global_recoil['Sr-88'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='r', label='Sr-88')
# ax.loglog(energy_group[:-1], global_recoil['Y-91'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='k', label='Y-91')

# ax.set_xlim(1e-6, 30)
# ax.set_ylim(1e+6, 1e+12)
# ax.legend()
# ax.grid()
# plt.show()

pkafig = SpectrumFig(title='PKA spectrum')
pkafig.add_plot(energy_group[:-1], global_recoil['H-1'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='g', label='H-1')
pkafig.add_plot(energy_group[:-1], global_recoil['Zr-92'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='y', label='Zr-92')
pkafig.add_plot(energy_group[:-1], global_recoil['Zr-93'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='k', label='Zr-93')
pkafig.add_plot(energy_group[:-1], global_recoil['Zr-95'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='b', label='Zr-95')
pkafig.add_plot(energy_group[:-1], global_recoil['Zr-97'].recoil_pka_spectrum * atoms_per_mole * density / atomic_mass, color='m', label='Zr-97')
pkafig.s_draw()
