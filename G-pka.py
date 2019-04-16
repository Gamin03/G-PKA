#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Main for G-pka.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-11 (last update: 2019-04-16)
@note   Rewritten the old G-PKA.
        - Add class definitions
        - Add the input file.
        - Change the pka xs matrix file into SPECTER-PKA format.
"""
from nuclide import Element, Nuclide, NuclideRecoil
from input import Input
from read_pka_file import *
from utility_pka import *
from utility_fig import *
from utility_basic import *
from utility_output import *

# Code version control
_major_version = 0
_minor_version = 1
_bugfix_version = 1
print_code_title_version(_major_version, _minor_version, _bugfix_version)

# Global variables
atoms_per_mole = 0.6022e+24

global_recoil = {}
global_element = {}

# --- Initialization ----------------------------------------------------------
if len(sys.argv) == 1:
    input_file = 'input.json'
    output_file = 'output.txt'
elif len(sys.argv) == 2:
    input_file = str(sys.argv[1])
    output_file = 'output.txt'
elif len(sys.argv) == 3:
    input_file = str(sys.argv[1])
    output_file = str(sys.argv[2])
else:
    raise ValueError('Input parameters ERROR!')

inp = Input(input_file)
# output = open(output_file, 'w')
output = sys.stdout
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

    # Deal each nuclide 处理每个原始靶核
    for k, recoil_nuc in enumerate(nuc.recoil_nuclides_particles_info):

        # Deal each reaction channel 处理每个反应道
        # --- Set up the NuclideRecoil for this nuc
        za_recoil = get_daughter_nuclides_particles(recoil_nuc[0], recoil_nuc[1], nuc.Z, nuc.A)
        nuc_recoil = NuclideRecoil(za_recoil[0], za_recoil[1], recoil_nuc[1])
        nuc_recoil.title = recoil_nuc[0]
        # --- Load the recoil matrix
        nuc_recoil.load_recoil_matrix(recoil_nuc[2])
        # --- Calculate recoil pka spectra and save it
        nuc_recoil.compute_recoil_pka_spectra(nuc.recoil_flux_pka)

        # Calculate dpa values
        if inp.do_damage:
            nuc_recoil.damage_function_coeffs, nuc_recoil.estimate_ed = \
                get_damage_coeffs_array(nuc.recoil_energy_group_struc, (nuc_recoil.Z, nuc_recoil.A), (nuc.Z, nuc.A))
            nuc_recoil.damage_cross_section = nuc_recoil.damage_function_coeffs * nuc_recoil.pka_spectrum

            nuc_recoil.damage_dpa = nuc_recoil.damage_cross_section * 0.8 / (2. * nuc_recoil.estimate_ed)

        nuc.append_recoil_nuclide(nuc_recoil)

        # Accumulate the NuclideRecoil for this nuc
        recoil_pka_spectrum_total += nuc_recoil.pka_spectrum

    # Save the total recoil pka spectrum for this nuc
    nuc.recoil_pka_spectrum = recoil_pka_spectrum_total

# --- Calculate total results ----------------------------------------------------------
# Check each recoil nuclides of every initial nuclide, save them into global
# nuclides and elements dictionary
# 计算每个初始靶核的每个反冲核，并保存在全局的核素和元素字典内
# ****** Add the nuclide ratio here ******
# ****** 注意，在此处引入核素比值 ******
print("\n\n>>> START CALCULATE TOTAL PKA RESULTS ...", file=output)
for nuc in nuclides:
    ratio = nuc.ratio

    for recoil_nuc in nuc.recoil_nuclides:
        # 对于每个反冲核，检查其核素或者元素在全局字典内是否存在
        # 若存在，则将当前pka谱加上；若不存在，则创建核素或元素，保存数据，并添加到全局字典内
        # Check if the nuclide exist in the dict, if not, append it
        if recoil_nuc.name in global_recoil:
            if ((recoil_nuc.name != 'He-4') and (800 <= recoil_nuc.mtd <= 849)) or \
               ((recoil_nuc.name != 'H-1')  and (600 <= recoil_nuc.mtd <= 649)):
                continue
            global_recoil[recoil_nuc.name].add_recoil_pka_spectrum(recoil_nuc.pka_spectrum * ratio)
            if inp.do_damage:
                global_recoil[recoil_nuc.name].add_damage_values(recoil_nuc, ratio)
        else:
            # Increase new nuclide into global nuclide dictionary
            nuclide = Nuclide((recoil_nuc.Z, recoil_nuc.A))
            nuclide.recoil_energy_group_struc = nuc.recoil_energy_group_struc
            nuclide.recoil_pka_spectrum = recoil_nuc.pka_spectrum * ratio
            if inp.do_damage:
                nuclide.copy_recoil_damage(recoil_nuc, ratio)
            global_recoil[recoil_nuc.name] = nuclide

        # Check if the element exist in the dict, if not, append it
        if recoil_nuc.element in global_element:
            if ((recoil_nuc.name != 'He-4') and (800 <= recoil_nuc.mtd <= 849)) or \
               ((recoil_nuc.name != 'H-1')  and (600 <= recoil_nuc.mtd <= 649)):
                continue
            global_element[recoil_nuc.element].add_recoil_pka_spectrum(recoil_nuc.pka_spectrum * ratio)
            if inp.do_damage:
                global_element[recoil_nuc.element].add_damage_values(recoil_nuc, ratio)
        else:
            # Increase new element into global element dictionary
            element = Element(recoil_nuc.Z)
            element.num_recoil_pka_energy_group = nuc.num_recoil_energy_group_struc
            element.copy_recoil_pka_energy_group(nuc.recoil_energy_group_struc)
            element.recoil_pka_spectrum = recoil_nuc.pka_spectrum * ratio
            if inp.do_damage:
                element.copy_recoil_damage(recoil_nuc, ratio)
            global_element[element.name] = element

# Calculate average pka and dpa values
# 计算平均的pka和dpa值
for key in global_recoil:
    nuclide = global_recoil[key]
    # in eV
    if sum(nuclide.recoil_pka_spectrum) > 0.:
        nuclide.average_pka_energy = sum(nuclide.recoil_pka_spectrum * 0.5 * 1.E+6 *
                        (nuclide.recoil_energy_group_struc[:-1] + nuclide.recoil_energy_group_struc[1:])) / \
                        sum(nuclide.recoil_pka_spectrum)

    nuclide.average_displacement_energy = sum(nuclide.damage_cross_section) * 1.E+6

for key in global_element:
    element = global_element[key]
    # in eV
    element.average_pka_energy = sum(element.recoil_pka_spectrum * 0.5 * 1.E+6 *
                            (element.recoil_pka_energy_group[:-1] + element.recoil_pka_energy_group[1:])) / \
                            sum(element.recoil_pka_spectrum)
    element.average_displacement_energy = sum(element.damage_cross_section) * 1.E+6

# ------ Sum results ------
total_pka_spectrum = None
for key in global_element:
    ele = global_element[key]
    # ele.calculate_average_pka_energy()
    if total_pka_spectrum is None:
        total_pka_spectrum = ele.recoil_pka_spectrum
    else:
        total_pka_spectrum += ele.recoil_pka_spectrum


# --- 输出结果 ----------------------------------------------------------------
# 将结果写入Excel表
if inp.do_write_each_nuclides:
    for nuc in nuclides:
        write_each_recoil_pka_into_xls(nuc)
write_total_nuclides_into_xls(global_recoil)
write_total_elements_into_xls(global_element)

print("\n\n>>> OUTPUT TOTAL PKA RESULTS ...", file=output)

output.close()

# 结果绘图
if inp.plot_figure:
    density = inp.density
    atomic_mass = inp.atomic_mass
    plot_global_element_figure(global_element, energy_group[:-1], [], density, atomic_mass, 'Element.png')
    plot_global_element_figure(global_recoil, energy_group[:-1], [], density, atomic_mass, 'Nuclide.png')
