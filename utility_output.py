#!/usr/bin/env python
"""
Utility output functions for the code.

@author Jimin Ma  <majm03@yeah.net>
@time   2018-09-20

"""

import sys
from xlwt import *

_head = ["Group Num", "Recoil energy (low)", "Recoil energy (high)", "PKAs", "PKAs norm_sum",
         "disp cross section", "NRT_dpa"]
# 创建配置，样式
_al = Alignment()
_al.horz = Alignment.HORZ_CENTER
_al.vert = Alignment.VERT_CENTER
_borders = Borders()
_borders.bottom = Borders.THICK
_style = XFStyle()
_style.alignment = _al
_style.borders = _borders
_style_num = XFStyle()
_style_num.num_format_str = '0.0000E+00'


# 将最终结果写入Excel文件，每个核素一个excel文件
def write_each_recoil_pka_into_xls(nuc):
    book = Workbook()
    # 对每个反冲核创建一个sheet
    for recoil in nuc.recoil_nuclides:
        sum_pka = sum(recoil.pka_spectrum)
        if sum_pka > 0:
            sheet = book.add_sheet('{},{}'.format(recoil.name, recoil.title[:-7]))
            # 写入标题行
            for i, text in enumerate(_head):
                sheet.row(0).write(i, text, style=_style)
            # 写入数据
            row = recoil.pka_spectrum.shape[0]
            row_values = 1
            for i in range(row):
                if recoil.pka_spectrum[i] > 0:
                    sheet.write(row_values, 0, i)
                    sheet.write(row_values, 1, nuc.recoil_energy_group_struc[i], _style_num)
                    sheet.write(row_values, 2, nuc.recoil_energy_group_struc[i+1], _style_num)
                    sheet.write(row_values, 3, recoil.pka_spectrum[i], _style_num)
                    sheet.write(row_values, 4, recoil.pka_spectrum[i] / sum_pka, _style_num)
                    sheet.write(row_values, 5, recoil.damage_cross_section[i], _style_num)      # damage - dpa
                    sheet.write(row_values, 6, recoil.damage_dpa[i], _style_num)                # damage - dpa
                    row_values += 1
            for i in range(len(_head)):
                sheet.col(i).width = 4000
    # 保存文件
    book.save('PKA_{}.xls'.format(nuc.name))


def write_total_nuclides_into_xls(global_recoil):
    book = Workbook()
    for key in global_recoil:
        print('[{}]'.format(key))
        sheet = book.add_sheet('{}'.format(key))
        nuclide = global_recoil[key]
        sum_pka = sum(nuclide.recoil_pka_spectrum)
        # 写入标题行
        for i, text in enumerate(_head):
            sheet.row(0).write(i, text, style=_style)
        # 写入数据
        row = nuclide.recoil_pka_spectrum.shape[0]
        row_values = 1
        for i in range(row):
            if nuclide.recoil_pka_spectrum[i] > 0:
                sheet.write(row_values, 0, i)
                sheet.write(row_values, 1, nuclide.recoil_energy_group_struc[i], _style_num)
                sheet.write(row_values, 2, nuclide.recoil_energy_group_struc[i + 1], _style_num)
                sheet.write(row_values, 3, nuclide.recoil_pka_spectrum[i], _style_num)
                sheet.write(row_values, 4, nuclide.recoil_pka_spectrum[i] / sum_pka, _style_num)
                sheet.write(row_values, 5, nuclide.damage_cross_section[i], _style_num)     # damage - dpa
                sheet.write(row_values, 6, nuclide.damage_dpa[i], _style_num)               # damage - dpa
                row_values += 1
        for i in range(len(_head)):
            sheet.col(i).width = 4000
        # 写入汇总行
        sheet.write(row_values, 1, 'Average PKA energy')
        sheet.write(row_values, 3, '{0:10.4E} (eV)'.format(nuclide.average_pka_energy))
        row_values += 1
        sheet.write(row_values, 1, 'Displacement energy')
        sheet.write(row_values, 3, '{0:10.4E} eV/s'.format(nuclide.average_displacement_energy))
        row_values += 1
        sheet.write(row_values, 1, 'Equivalent NRT dpa')
        # **** 注意这里的40. *** 需修改
        sheet.write(row_values, 3, '{0:10.4E} dpa/s'.format(nuclide.average_displacement_energy * 0.8 / (2. * 40.)))
    # 保持文件
    book.save('Total_PKAs_nuclides.xls')


def write_total_elements_into_xls(global_element):
    book = Workbook()
    for key in global_element:
        sheet = book.add_sheet('{}'.format(key))
        element = global_element[key]
        sum_pka = sum(element.recoil_pka_spectrum)
        # 写入标题行
        for i, text in enumerate(_head):
            sheet.row(0).write(i, text, style=_style)
        # 写入数据
        row = element.recoil_pka_spectrum.shape[0]
        row_values = 1
        for i in range(row):
            if element.recoil_pka_spectrum[i] > 0:
                sheet.write(row_values, 0, i)
                sheet.write(row_values, 1, element.recoil_pka_energy_group[i], _style_num)
                sheet.write(row_values, 2, element.recoil_pka_energy_group[i + 1], _style_num)
                sheet.write(row_values, 3, element.recoil_pka_spectrum[i], _style_num)
                sheet.write(row_values, 4, element.recoil_pka_spectrum[i] / sum_pka, _style_num)
                sheet.write(row_values, 5, element.damage_cross_section[i], _style_num)         # damage - dpa
                sheet.write(row_values, 6, element.damage_dpa[i], _style_num)                   # damage - dpa
                row_values += 1
        for i in range(len(_head)):
            sheet.col(i).width = 4000
        # 写入汇总行
        sheet.write(row_values, 1, 'Average PKA energy')
        sheet.write(row_values, 3, '{0:10.4E} (eV)'.format(element.average_pka_energy))
        row_values += 1
        sheet.write(row_values, 1, 'Displacement energy')
        sheet.write(row_values, 3, '{0:10.4E} eV/s'.format(element.average_displacement_energy))
        row_values += 1
        sheet.write(row_values, 1, 'Equivalent NRT dpa')
        # **** 注意这里的40. *** 需修改
        sheet.write(row_values, 3, '{0:10.4E} dpa/s'.format(element.average_displacement_energy * 0.8 / (2. * 40.)))
    # 保存文件
    book.save('Total_PKAs_elements.xls')
