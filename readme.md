# G-PKA

> 作者：马纪敏 [majm03@yeah.net](majm03@yeah.net)

## 概述

该程序用来计算核素、元素或化合物在中子或质子通量辐照下的PKA能谱、dpa损伤截面。

## 运行方式

### 必需文件

运行G-PKA必需的文件有：

- G-PKA可执行程序 G-PKA.exe；

- 输入文件，默认文件名为 input.json，格式见后说明；

- 能谱数据文件，

- 核素数据文件

运行方式为，在命令行内输入如下命令
```bash
G-pka.exe [input.json] [output.txt]
```

### 输入文件

G-PKA的输入文件采用 json 格式。
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