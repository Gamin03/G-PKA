#!/usr/bin/env python
"""
Utility figure functions for the code.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-17

"""

import matplotlib.pyplot as plt
import sys

_avogadro = 6.022141930E+23


class SpectrumFig:
    def __init__(self, title):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        self.args = []
        self.kwargs = []

    def add_plot(self, *args, **kwargs):
        self.args.append(args)
        self.kwargs.append(kwargs)

    def s_draw(self, figure_name):
        for args, kwargs in zip(self.args, self.kwargs):
            self.ax.loglog(*args, **kwargs)

        self.ax.set_xlim(1e-6, 30)
        self.ax.set_ylim(1e+6, 1e+14)
        self.ax.legend()
        self.ax.grid(color='grey', which='major', linestyle=':')
        self.ax.set_xlabel('PKA energy ($MeV$)')
        self.ax.set_ylabel(r'PKAs $s^{-1} cm^{-3}$')
        # plt.show()
        plt.savefig(figure_name)


# 绘制所有的核素和元素图
def plot_global_element_figure(global_dicts, e_group, pick_list, density, atomic_mass, figure_title):
    linestyles = ['-', '-.', ':']
    colors = ["red", "yellow", "blue", "lightgreen", "black"]
    num = len(pick_list)
    if num == 0:
        num = len(global_dicts)
    pkafig = SpectrumFig(title='PKA spectrum')
    # 若pick_list为空，则绘制所有字典内的数据
    if len(pick_list) == 0:
        if len(global_dicts) < 11:
            idx = 0
            for key in global_dicts:
                # plot figure
                pkafig.add_plot(e_group, global_dicts[key].recoil_pka_spectrum * _avogadro * density / atomic_mass,
                                linestyle=linestyles[idx % 3], color=colors[idx % 5], label='{}'.format(key))
                idx += 1
        else:
            # 若字典多于10个，则挑选最大的10个进行绘图
            global_lists_sorted = sorted(dict_to_list(global_dicts), key=lambda x: sum(x[1].recoil_pka_spectrum),
                                         reverse=True)
            idx = 0
            for i in range(10):
                pkafig.add_plot(e_group, global_lists_sorted[i][1].recoil_pka_spectrum * _avogadro * density / \
                                atomic_mass, linestyle=linestyles[idx % 3], color=colors[idx % 5],
                                label='{}'.format(global_lists_sorted[i][0]))
                idx += 1
                if idx == 10:
                    break

    else:
        # Check if pick_list are in global_dicts key
        for i in range(len(pick_list)):
            if not global_dicts.has_key(pick_list[i]):
                print("ERROR in pick_list: {} is not in global_nuclide dict!".format(pick_list[i]))
                sys.exit()
        # plot figure
        for i in range(len(pick_list)):
            key = pick_list[i]
            pkafig.add_plot(e_group, global_dicts[key].recoil_pka_spectrum * _avogadro * density / atomic_mass,
                            linestyle=linestyles[i % 3], color=colors[i % 5], label='{}'.format(key))

    pkafig.s_draw(figure_title)


# 将字典转化为列表
def dict_to_list(dic: dict):
    keys = dic.keys()
    vals = dic.values()
    lst = [(key, val) for key, val in zip(keys, vals)]
    return lst
