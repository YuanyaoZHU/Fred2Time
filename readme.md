# 本程序来源

- 本程序是从YingyiLiu的Fred2Timd工程fork过来的。

- 该程序由本文对求解器进行了重写，采用的是lungkuta4阶方法。

- 改程序存在bug，即release版本的.exe程序在计算时会出现系泊无法收敛的问题，在debug模式下却可行，初步判断可能时编译
器的问题。

- **Inertia Property.xlsx**文件可以用于计算浮体的惯性张量矩阵。

- **Wamit2Hams.m**文件是用来将Wamit水动力文件转换成Hams水动力文件格式的matlab程序。