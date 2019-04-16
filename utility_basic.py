#!/usr/bin/env python
"""
Utility functions for the code.

@author Jimin Ma  <majm03@yeah.net>
@time   2018-09-19

"""

import sys


# 程序标题及版本号输出
def print_code_title_version(major_version, minor_version, bugfix_version, file=sys.stdout):
    print("--------------------------------------------------", file=file)
    print("G-PKA v{}.{}.{}".format(major_version, minor_version, bugfix_version), file=file)
    print("", file=file)
    print("    Author   : Jimin Ma <majm03@yeah.net>", file=file)
    print("    Copyright: INPC, CAEP @ Sichuan, China", file=file)
    print("--------------------------------------------------\n", file=file)
