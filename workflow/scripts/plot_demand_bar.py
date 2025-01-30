#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 09:32:38 2023

@author: fabian
"""

import os
import pypsa
import pypsa.statistics
import seaborn as sns
from common import import_network, mock_snakemake
import matplotlib.pyplot as plt

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_demand_bar", ext="png", clusters=90, run="baseline"
    )


sns.set_theme(**snakemake.params["theme"])

n = import_network(snakemake.input.network)
s = n.statistics

demand = s.withdrawal(groupby=pypsa.statistics.get_carrier_and_bus_carrier).Load


rename_bus_carrier = {
    "solid biomass for industry": "biomass",
    "gas for industry": "gas",
    "low voltage": "electricity",
    "heat": "low-t heat",
    "Li ion": "electricity",
    "oil": "liquid\nhydrocarbons",
}
rename_carrier = {
    "Heat": "Residential and Services Heat",
    "Electricity": "Residential and Services Elec.",
}

demand = (
    demand.rename(rename_bus_carrier, level=1)
    .rename(rename_carrier, level=0)
    .rename(str.title, level=1)
    .rename(n.carriers.nice_name, level=0)
)
carrier_order = demand.loc[
    :, demand.groupby(level=1).sum().sort_values().index
].index.get_level_values(0)


df = demand.unstack().div(1e9)
df = df[df.sum().sort_values().index].loc[carrier_order]

color = n.carriers.set_index("nice_name").rename(rename_carrier).color.to_dict()


fig, ax = plt.subplots(figsize=(5, 4))
df.T.plot.bar(stacked=True, color=color, xlabel="", ylabel="Energy demand [PWh]", ax=ax)
pad = 0.01
for i, val in enumerate(df.sum()):
    ax.text(i, val + pad, f"{val:.2f}", ha="center", va="bottom", fontsize=8, color="k")

ax.legend(title="", frameon=False, bbox_to_anchor=(0, 1.2), loc="upper left")
sns.despine()

fig.savefig(
    snakemake.output.figure,
    bbox_inches="tight",
)
df.fillna(0).round(3).to_csv(snakemake.output.table)
