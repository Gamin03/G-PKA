#!/usr/bin/env python
"""
Class definition for input file.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-11

"""

import os
import json
import numpy as np


class Input:
    def __init__(self, name):
        self.infile_name = name
        self.infile_flux_name = None
        self.number_pka_files = None
        self.pka_files = None
        self.flux_rescale_value = 1.0
        self.assumed_ed = 40.0

        self.flux_unit = None   # 0 is n s^-1^, 1 is n s^-1^ MeV^-1^
        self.num_flux_energy_group = 0
        self.flux_energy_group = None
        self.flux_spectrum = None

        self.do_gamma_estimate = False

    def read_infile(self, file_object):
        with open(self.infile_name) as f:
            print(">>> START READ INPUT FILE [{}]".format(self.infile_name), file=file_object)
            data = json.load(f)
            self.infile_flux_name = data['flux_filename']
            self.number_pka_files = data["number_pka_files"]
            self.pka_files = data["columns"]
            self.flux_rescale_value = data["flux_rescale_value"]
            self.assumed_ed = data["assumed_ed"]
            self.do_gamma_estimate = data["do_gamma_estimate"]
            print("--- FINISH READING INPUT FILE [{}]\n".format(self.infile_name), file=file_object)

    def read_flux(self, file_object):
        if self.infile_flux_name is None:
            print("No input flux file!")
        with open(self.infile_flux_name) as flux_file:
            print(">>> START READ FLUX INPUT FILE [{}]".format(self.infile_flux_name), file=file_object)
            flux_file.readline()
            line = flux_file.readline()  # 跳过标题行和第二行
            line = line.strip().split()
            if int(line[2]) == 2:
                self.flux_unit = 'n s^{-1}'
            else:
                self.flux_unit = 'n s^{-1} MeV^{-1}'
            line = flux_file.readline()
            line = line.strip().split()
            group = self.num_flux_energy_group = int(line[0])
            self.flux_spectrum = np.zeros(group)
            self.flux_energy_group = np.zeros(group + 1)
            for i in range(group + 1):
                line = flux_file.readline()
                self.flux_energy_group[i] = float(line.strip())
            for i in range(group):
                line = flux_file.readline()
                self.flux_spectrum[i] = float(line.strip())

        # ------ Rescale the flux spectrum ------
        if self.flux_unit == 'n s^{-1} MeV^{-1}':
            # Change unit into 'n s^{-1}'
            self.flux_spectrum *= (self.flux_energy_group[1:] - self.flux_energy_group[:-1])
        # Make sure the normalization
        total_flux = self.flux_spectrum.sum()
        self.flux_spectrum /= total_flux
        # Rescale
        self.flux_spectrum *= self.flux_rescale_value
        # Make sure the flux spectrum unit is 'n s^{-1} MeV^{-1}'
        self.flux_spectrum /= (self.flux_energy_group[1:] - self.flux_energy_group[:-1])
        self.flux_unit = 'n s^{-1} MeV^{-1}'

        # ------ Output flux information ------
        print("\tFlux group number  : {}".format(self.num_flux_energy_group), file=file_object)
        print("\tTotal flux         : {0:.3e} {1}".format(self.flux_rescale_value, self.flux_unit), file=file_object)


if __name__ == '__main__':
    ip = Input('input.json')
    ip.read_infile()

    print(ip.pka_files[0]['parent'])

