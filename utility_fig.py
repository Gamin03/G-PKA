#!/usr/bin/env python
"""
Utility figure functions for the code.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-17

"""

import matplotlib.pyplot as plt

class SpectrumFig:
    def __init__(self, title):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

        self.args = []
        self.kwargs = []

    def add_plot(self, *args, **kwargs):
        self.args.append(args)
        self.kwargs.append(kwargs)

    def s_draw(self):
        for args, kwargs in zip(self.args, self.kwargs):
            self.ax.loglog(*args, **kwargs)

        self.ax.set_xlim(1e-6, 30)
        self.ax.set_ylim(1e+6, 1e+12)
        self.ax.legend()
        self.ax.grid(color='grey', which='major', linestyle=':')
        self.ax.set_xlabel('PKA energy ($MeV$)')
        self.ax.set_ylabel(r'PKAs $s^{-1} cm^{-3}$')
        plt.show()
