#!/usr/bin/env python
"""
Class definition for input file.

@author Jimin Ma  <majm03@foxmail.com>
@time   2018-09-11

"""

import numpy as np
import scipy.sparse as sp


# Read 1st part in pka file (SPECTER-PKA pka-file format)
# Get the pka incident and recoil energy structure
def read_pka_file_energy_group_struc(file_object):
    # Check the title line
    line = file_object.readline()

    title = line[:30].strip()
    mtd = int(line[30:35])
    line_right = line[36:].strip().split()
    num_pka_incident_egs = int(line_right[1]) - 1
    num_pka_recoil_egs = num_pka_incident_egs
    num_pka_points = int(line_right[2])
    assert (num_pka_incident_egs == num_pka_points), "PKA points must be same as incident energy group number!"

    # Read the pka spectrum
    i_label = 0
    pka_e_array = np.zeros((num_pka_incident_egs + 1) * 2)
    for i in range((num_pka_incident_egs + 1) * 2 // 6):  # floor division
        line = file_object.readline()
        line = line.strip().split()
        for j in range(6):
            pka_e_array[i_label] = float(line[j])
            i_label += 1

    # Read the last spectrum line (if exist)
    tmp = (num_pka_incident_egs + 1) * 2 % 6
    if tmp > 0:
        line = file_object.readline()
        line = line.strip().split()
        for j in range(tmp):
            pka_e_array[i_label] = float(line[j])
            i_label += 1

    # Read the last line - ONLY for the 1st part
    # Read the sparse matrix
    ng = num_pka_incident_egs
    rows = []
    columns = []
    values = []
    while True:
        x = file_object.tell()
        line = file_object.readline()
        if 'matrix' in line or 'section' in line:
            file_object.seek(x)  # Back to last title line
            break
        line = line.strip().split()
        rows.append(int(line[0]) - 1)
        columns.append(int(line[1]) - 1)
        values.append(float(line[2]))
    # Save matrix
    M = sp.coo_matrix((np.array(values), (np.array(rows), np.array(columns))), shape=[ng, ng])

    return title, mtd, num_pka_incident_egs, pka_e_array[:(num_pka_incident_egs+1)], M


# Read each recoil matrix section
def read_pka_file_each_matrix(file_object):
    # Read the title line
    x0 = file_object.tell()
    line = file_object.readline()
    title = line[:30].strip()
    mtd = int(line[30:35])
    line_right = line[36:].strip().split()
    ng = line_right[2]

    if 'matrix' in title:
        # Read the sparse matrix
        rows = []
        columns = []
        values = []
        while True:
            x = file_object.tell()
            line = file_object.readline()
            if 'matrix' in line or 'section' in line:
                file_object.seek(x)  # Back to last title line
                break
            line = line.strip().split()
            rows.append(int(line[0]) - 1)
            columns.append(int(line[1]) - 1)
            values.append(float(line[2]))
        # Save matrix
        M = sp.coo_matrix((np.array(values), (np.array(rows), np.array(columns))), shape=[ng, ng])
    else:
        file_object.seek(x0)
        M = None  # if matrix don't exist

    return title, mtd, M


# Read the (n, gamma) cross section
def read_pka_file_ng_xs(file_object):
    line = file_object.readline()
    title = line[:30].strip()
    mtd = int(line[30:35])
    line_right = line[36:].strip().split()
    ng = int(line_right[1]) - 1

    ng_xs = np.zeros(ng)
    while True:
        line = file_object.readline()
        if line == "":
            break
        line = line.strip().split()
        ng_xs[int(line[0]) - 1] = float(line[1])

    return title, mtd, ng_xs


if __name__ == '__main__':
    with open('./test/example_data/F019s-p.asc') as fp:
        title, mtd, pi, pe_array, A = read_pka_file_energy_group_struc(fp)
        print("{0:30s} {1:4d}".format(title, mtd))
        while True:
            title, mtd, A = read_pka_file_each_matrix(fp)
            if A is None:
                break
            else:
                print("{0:30s} {1:4d}".format(title, mtd))
        title, mtd, ng_xs = read_pka_file_ng_xs(fp)
        print("{0:30s} {1:4d}".format(title, mtd))


