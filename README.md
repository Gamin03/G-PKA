#  G-PKA

Author: Jimin MA 

E-mail: *majm03@yeah.net*

Time: 2019-04-16

Version: v0.1.2


##  Introduction

The code is used to calculate the PKA spectrum and dpa damage cross section of nuclide, element or compound under neutron or proton irradiation.

Written by Python3. Required packages are:

numpy>=1.15.0, scipy>=1.1.0, xlwt>=1.3.0, json>=2.6.0, matplotlib>=3.0.0

Please cited this paper if needed:

**Jimin Ma, Hongwen Huang. PKA spectra and irradiation damage calculations. 8th Reactor Physics and Nuclear Material Conference. Shenzhen, China, 2017**

##  Run method

### Requirement

- input file, default name is 'input.json'. Detail of the input format see next part.
- spectrum data file.
- nuclide data file. Use same file as SPECTRA-PKA, which can be download in https://fispact.ukaea.uk/nuclear-data/downloads/ .

  Run the code in terminal as follows: 
```python
python G-pka.py [input.json] [outout.json]
```

### Input file

The G-PKA input file uses *json* format.

Main parameters are:

| Parameter name     | value meaning                                      |
| ------------------ | -------------------------------------------------- |
| flux_filename      | flux file name. The format is same as SPECTER-PKA. |
| number_pka_files   | pka data file counts.                              |
| columns            | pka file and nuclide contents.                     |
| flux_rescale_value | flux factor.                                       |
| assumed_ed         | Assumed Ed value in dpa calculation.               |
| do_gamma_estimate  | Does gamma dose be estimated in pka calculation?   |
| plot_figure        | Plot figure option。                               |

Second level parameters in `columns` are:

| parameter name       | value meaning                                 |
| -------------------- | --------------------------------------------- |
| pka_filename         | pka data file name (include directory).       |
| pka_ratios           | pka nuclide cpntents in total nuclides.       |
| parent               | current pka nuclide name.                     |
| ngamma_parent_mass   | parent nuclide mass ( used in gamma estimate) |
| ngamma_daughter_mass | daughter mass ( used in gamma estimate)       |

### Result file

The detail pka and dpa values of nuclides are given in *excel* file in G-PKA calculation. Result file names for nuclides and elements are *Total_PKAs_nuclides.xls* and *Total_PKAs_elements.xls*, respectively.

The pka spectrum of nuclide and elements (the first larges 10) are plotted in figures.

### Example

The *input.json* file in the root directory is an example which calculate PKA for natural Zr under neutron irradiation in PWR spectrum. The example is same as that in SEPCTRA-PKA. Therequired data and result files are in *test* dirctory. 



Thanks! Hope the code helpful for you!

