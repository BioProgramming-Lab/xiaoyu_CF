import xiaoyu_CF as flow
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

my_expr = flow.xiaoyu_Expr("26Mar2025/HuaiDan-New Experiment-20250326-1722")
my_expr.add_conc_condition(1E-5, dilu_dir="down")

#de_outliers = flow.assist_gate(my_expr.expr, if_plot = True, low = 0.5E9, high = 1.75E9)
cell_type = flow.auto_gate_FSC_SSC(my_expr.expr, if_plot = True)
de_doublets = flow.auto_gate_FSC_A_H(cell_type, if_plot = True, keep = 0.7)

af_op = flow.cytoflow.AutofluorescenceOp()
af_op.channels = ["mCherry-A", "TagBFP-A", "549/15-488 nm citrine-A", "Alexa 647-A"]
af_op.blank_file = "21Apr2024/blank.fcs"
af_op.estimate(de_doublets)

ex_af = af_op.apply(de_doublets)

bl_op = flow.cytoflow.BleedthroughLinearOp()
bl_op.controls = {"TagBFP-A": "21Apr2024/testGroup/B1.fcs",
                  "mCherry-A": "21Apr2024/testGroup/B2.fcs",
                  "549/15-488 nm citrine-A": "21Apr2024/testGroup/B3.fcs",
                  "Alexa 647-A": "21Apr2024/testGroup/B7.fcs"}
bl_op.estimate(ex_af)

ex_bl = bl_op.apply(ex_af)

x_class = []
possitive_rate = []
for data_set in ['A', 'B', 'C']:
    mix = flow.subset_by_char(de_doublets, data_set)

    EpCAM = flow.assist_gate(mix, if_plot = True, channel = "mCherry-A", low = 5E6, high = 1E9)
    EGFR = flow.assist_gate(mix, if_plot = True, channel = "TagBFP-A", low = 4E6, high = 1E9)
    CoExpression = flow.assist_gate(mix, if_plot = True, channel = "549/15-488 nm citrine-A", low = 4E6, high = 1E9)

    EGFR_hist = EGFR.data.groupby('sample')['Alexa 647-A'].apply(list).tolist()[0]
    EpCAM_hist = EpCAM.data.groupby('sample')['Alexa 647-A'].apply(list).tolist()[0]
    CoExpression_hist = CoExpression.data.groupby('sample')['Alexa 647-A'].apply(list).tolist()[0]

    fig, ax = plt.subplots(figsize=(3, 2))
    sns.kdeplot(EGFR_hist, ax=ax, label="EGFR", fill=True, alpha=0.5)
    sns.kdeplot(EpCAM_hist, ax=ax, label="EpCAM", fill=True, alpha=0.5)
    sns.kdeplot(CoExpression_hist, ax=ax, label="CoExpression", fill=True, alpha=0.5)
    # sns.kdeplot(nega_hist, ax=ax, label="NegativeCtl", fill=True, alpha=0.5)
    ax.axvline(x=3E8, color='r', linestyle='--', label='Cut-off')
    plt.xlabel("Flourescence Intensity")
    plt.xlim(-1e8, 1.2e9)
    plt.xticks([0, 3e8, 6e8, 9e8, 1.2e9])
    plt.yticks([])
    ax.legend()
    plt.savefig(f"26Mar2025/mix_discrimination_lin_{data_set}.svg")
    #plt.show()

    EGFR_positive = flow.assist_gate(EGFR, if_plot = True, channel = "Alexa 647-A", low = 3E8, high = 1E10)
    EpCAM_positive = flow.assist_gate(EpCAM, if_plot = True, channel = "Alexa 647-A", low = 3E8, high = 1E10)
    CoExpression_positive = flow.assist_gate(CoExpression, if_plot = True, channel = "Alexa 647-A", low = 3E8, high = 1E10)

    print(data_set, ":\t", len(EGFR_positive)/len(EGFR), len(EpCAM_positive)/len(EpCAM), len(CoExpression_positive)/len(CoExpression))
    x_class.append("CHO/EGFR")
    possitive_rate.append(len(EGFR_positive)/len(EGFR))
    x_class.append("CHO/EpCAM")
    possitive_rate.append(len(EpCAM_positive)/len(EpCAM))
    x_class.append("CHO/EGFR/EpCAM")
    possitive_rate.append(len(CoExpression_positive)/len(CoExpression))
    '''
    plt.figure()
    plt.bar(["EGFR", "EpCAM", "CoExpression"], [len(EGFR_positive)/len(EGFR), len(EpCAM_positive)/len(EpCAM), len(CoExpression_positive)/len(CoExpression)])
    plt.title("Positive Rate")
    plt.savefig("21Apr2024/positive_rate.svg")
    plt.show()'''

# Create a DataFrame from the data
data = pd.DataFrame({'Cell Line': x_class, 'Positive Rate': possitive_rate})

# Calculate mean values for each group
data_means = data.groupby('Cell Line')['Positive Rate'].mean().reset_index()

ColorEpCAM = '#7277A1'
ColorEGFR = '#E77751'
ColorCoExpr = '#AB774C'

# Plotting the swarm plot with specified colors
plt.figure(figsize=(2.5, 2.5))  # Adjusting the width (5) and height (6) of the figure
palette = {'CHO/EGFR': ColorEGFR, 'CHO/EpCAM': ColorEpCAM, 'CHO/EGFR/EpCAM': ColorCoExpr}
ax = sns.swarmplot(x='Cell Line', y='Positive Rate', data=data, palette=palette)

# Adding horizontal lines for the mean values
x_unique = data['Cell Line'].unique()
for i, category in enumerate(x_unique):
    mean_value = data_means[data_means['Cell Line'] == category]['Positive Rate'].values[0]
    x_start = i - 0.2
    x_end = i + 0.1
    ax.plot([x_start, x_end], [mean_value, mean_value], color='black', linestyle='-', linewidth=1)
ax.set_ylim(-0.05, 1)  # Set the y-axis limits


plt.ylabel('Positive rate', fontsize=10)
plt.xticks(rotation=20, fontsize=10)  # Rotate x-axis labels by 45 degrees
plt.yticks(fontsize=10)

# Remove top and right spines
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)

plt.tight_layout()  # Adjust the layout to prevent label overlap
plt.savefig('26Mar2025/positive_rate.svg', format='svg')
plt.show()