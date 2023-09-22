#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 18:11:52 2023

@author: fabian
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from common import (
    import_network,
    mock_snakemake,
)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_energy_balance_area",
        ext="pdf",
        kind="oil",
        clusters=40,
        run="projected-price",
    )

sns.set_theme(**snakemake.params["theme"])

config = snakemake.config
kind = snakemake.wildcards.kind
kind_to_carrier = {"hydrogen": "H2", "carbon": "co2 stored", "electricity": "AC"}
carrier = kind_to_carrier.get(kind, kind)
label = config["labels"].get(kind, kind.title())


n = import_network(snakemake.input.network)
key = snakemake.config["labels"][n.meta["wildcards"]["run"]]
balance = n.statistics.energy_balance(aggregate_time=False)


norm = 1e6
unit = "TWh" if kind not in ["carbon", "co2"] else "Mt"

ds = balance.xs(carrier, level=2)
ds = ds.div(norm)
ds = ds.droplevel(0).dropna(how="all").fillna(0).groupby(level=0).sum()
ds = ds.sort_values(ds.columns[0], ascending=False)

fig, ax = plt.subplots(1, 1, figsize=snakemake.params.settings["figsize"])

colors = n.carriers.set_index("nice_name").color.to_dict()
ds.T.plot.bar(color=colors, ax=ax, stacked=True, alpha=0.8, lw=0, width=1, rot=0)


ax.axhline(0, color="k", lw=1)
ax.set_ylabel(f"{label} [{unit}]")
ax.grid(axis="y", alpha=0.5)
if snakemake.params.settings.get("title", True):
    ax.set_title(f"{label} Balance")

# Set 12 x-ticks evenly spaced along the x-axis
num = len(ds.columns)
ax.set_xticks(np.arange(0, num, num / 12))
ax.set_xticklabels(
    [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ],
    ha="left",
)
ax.set_xlabel("")

# legend control is nasty
by = ds.columns[0]
pindex = ds[ds.gt(0).any(axis=1)].sort_values(by, ascending=True).index
nindex = ds[ds.le(0).any(axis=1)].sort_values(by, ascending=False).index
h, l = ax.get_legend_handles_labels()
order = [*pindex, *[i for i in nindex if i not in pindex]]
legend = pd.Series(h, index=l).drop_duplicates().reindex(order)
ax.legend(
    legend.values,
    legend.index,
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
    ncol=1,
)

sns.despine()
fig.savefig(snakemake.output.figure, dpi=300, bbox_inches="tight")
