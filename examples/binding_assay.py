import marimo

__generated_with = "0.23.8"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md("""
    # Cell Binding Assay — anti-His
    """)
    return


@app.cell
def _():
    import sys
    sys.path.insert(0, "examples")

    import warnings
    warnings.filterwarnings("ignore")

    import matplotlib.pyplot as plt
    import xiaoyu_CF as flow

    data_dir = "test/HuaiDan-cell binding assay anti-His-20240417-2003"
    expr = flow.xiaoyu_Expr(data_dir)
    print(f"Events: {len(expr.expr):,}")

    fsc = flow._find_channel(expr.expr, "FSC", "-A")
    ssc = flow._find_channel(expr.expr, "SSC", "-A")
    fsc_h = flow._find_channel(expr.expr, "FSC", "-H")
    return expr, flow, fsc, ssc


@app.cell
def _(mo):
    mo.md("""
    ## 1. Cleanup — remove saturation
    """)
    return


@app.cell
def _(expr, flow):
    clean = flow.gate_saturation(expr.expr)
    print(f"Removed {len(expr.expr) - len(clean):,} events ({len(clean):,} remain)")
    return (clean,)


@app.cell
def _(mo):
    mo.md("""
    ## 2. Gate Cells
    """)
    return


@app.cell
def _(clean, flow, mo):
    cells, _fig = flow.auto_gate_FSC_SSC(clean, keep=0.6, return_fig=True)
    print(f"Cells: {len(cells):,} ({len(cells)/len(clean)*100:.0f}%)")
    mo.mpl.interactive(_fig)
    return (cells,)


@app.cell
def _(mo):
    mo.md("""
    ## 3. Gate Singlets
    """)
    return


@app.cell
def _(cells, flow, mo):
    single, _fig = flow.auto_gate_FSC_A_H(cells, keep=0.7, return_fig=True)
    print(f"Singlets: {len(single):,}")
    mo.mpl.interactive(_fig)
    return (single,)


@app.cell
def _(mo):
    mo.md("""
    ## 4. Histogram — APC
    """)
    return


@app.cell
def _(mo):
    mo.md("""
    ## 4. Manual Gate (draw polygon, then re-run next cell)
    """)
    return


@app.cell
def _(clean, flow, fsc, mo, ssc):
    _fig = flow.interactive_gate_preview(clean, xchannel=fsc, ychannel=ssc,
                                           xscale="linear", yscale="linear")
    mo.mpl.interactive(_fig)
    return


@app.cell
def _(clean, flow):
    # re-run this cell after drawing the polygon above
    hand_gated = flow.apply_drawn_gate(clean, name="hand_gate")
    print(f"Manual gate: {len(hand_gated):,} events")
    return


@app.cell
def _(mo, single):
    import flow_analysis as fa
    _fig = fa.plot_histogram(single, channel="APC (Allophycocyanin)-A", huefacet="sample")
    mo.mpl.interactive(_fig)
    return (fa,)


@app.cell
def _(mo):
    mo.md("""
    ## 5. Ridgeline — APC
    """)
    return


@app.cell
def _(fa, mo, single):
    _fig = fa.plot_ridgeline(single, channel="APC (Allophycocyanin)-A", huefacet="sample")
    mo.mpl.interactive(_fig)
    return


@app.cell
def _(mo):
    mo.md("""
    ## 6. 96-well Median
    """)
    return


@app.cell
def _(expr, mo, single):
    med = expr.median_96well(single, channel="APC (Allophycocyanin)-A")
    import pandas as pd
    _df = pd.DataFrame(med, index=list("ABCDEFGH"), columns=[str(c) for c in range(1, 13)])
    mo.ui.table(_df.round(0))
    return


if __name__ == "__main__":
    app.run()
