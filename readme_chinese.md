# G-PKA 

> 作者：马纪敏 [majm03@yeah.net](majm03@yeah.net)
>
> 时间：2019-04-16
>
> 版本：v0.1.1

## 概述

该程序用来计算核素、元素或化合物在中子或质子通量辐照下的PKA能谱、dpa损伤截面。

Python3程序，需要的包有：

numpy>=1.15.0, scipy>=1.1.0, xlwt>=1.3.0, json>=2.6.0, matplotlib>=3.0.0

若需要引用，请引用下文：

<u>马纪敏，黄洪文. 中子辐照下的PKA能谱及辐照损伤计算. 第八届反应堆物理与核材料学术讨论会. 深圳，2017 </u>



## 运行方式

### 必需文件

运行G-PKA必需的文件有：

- G-PKA相关源码；

- 输入文件，默认文件名为 input.json，格式见后说明；

- 能谱数据文件，

- 核素数据文件，与SPECTRA-PKA相同，可在https://fispact.ukaea.uk/nuclear-data/downloads/下载。

运行方式为，在终端内输入如下命令（括号内为可选的默认参数）
```bash
python G-pka.py [input.json] [output.txt]
```

### 输入文件

G-PKA的输入文件采用 json 格式。默认输入文件名为input.json。

参数为

| 参数            | 内容                      |
|-----------------|---------------------------|
| flux_filename   | 通量文件名，格式同SPECTER。 |
| number_pka_files| pka数据文件数目。          |
| columns         | pka文件及核素组成说明。     |
| flux_rescale_value | 通量因子。              |
| assumed_ed      | 假定的dpa中Ed值。          |
| do_gamma_estimate | pka计算是否进行gamma剂量估算。 |
| plot_figure     | 是否进行绘图。             |

其中 `columns`中的二级参数为

| 参数            | 内容                        |
|-----------------|-----------------------------|
| pka_filename    | pka数据文件名称（含目录）。     |
| pka_ratios      | pka核素占总核素份额。         |
| parent          | 当前pka核素名称。            |
| ngamma_parent_mass | 母核质量（gamma估算使用）。 |
| ngamma_daughter_mass | 子核质量（gamma估算使用）。|

### 结果文件

G-PKA 计算的详细核素pka和dpa值通过excel文档给出。
核素为Total_PKAs_nuclides.xls，元素为Total_PKAs_elements.xls。

元素和核素（最大的10个）的pka能谱通过绘图给出。

### 算例

根目录下的input.json为计算Zr的PKA算例。直接运行即可。算例与SPECTRA-PKA的算例一致。所需的数据和结果文件在test文件夹内。

### 

