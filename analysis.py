import xiaoyu_CF as flow
import matplotlib.pyplot as plt

my_expr = flow.xiaoyu_Expr("21Mar2024/HuaiDan-cell binding assay-20240321-1952")
my_expr.add_conc_condition(1E-5, dilu_dir="down")

de_outliers = flow.gate_outliers(my_expr.expr)
cell_type = flow.auto_gate_FSC_SSC(de_outliers, if_plot = True)
de_doublets = flow.auto_gate_FSC_A_H(cell_type, if_plot = True, keep = 0.7)

flow.channel_histogram(flow.subset_by_num(de_doublets, '12'), huefacet="conc")

my_expr.median_96well(de_doublets)

plt.show()