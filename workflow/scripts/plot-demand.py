#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 09:32:38 2023

@author: fabian
"""

import pypsa
import seaborn as sns
from common import import_network
import matplotlib.pyplot as plt

sns.set_theme(style="white", context="paper")

n = import_network(
    "/home/fabian/papers/co2-network/results/baseline/networks/90_nodes.nc"
)
s = n.statistics


demand = s.withdrawal(groupby=s.groupers.get_carrier_and_bus_carrier).Load


rename_bus_carrier = {
    "solid biomass for industry": "biomass",
    "gas for industry": "gas",
    "low voltage": "electricity",
    "services rural heat": "Low-Temperature\nHeat",
    "urban central heat": "Low-Temperature\nHeat",
    "residential rural heat": "Low-Temperature\nHeat",
    "residential urban decentral heat": "Low-Temperature\nHeat",
    "services urban decentral heat": "Low-Temperature\nHeat",
    "Li ion": "electricity",
}


demand = (
    demand.rename(rename_bus_carrier, level=1)
    .rename(str.title, level=1)
    .rename(n.carriers.nice_name, level=0)
)
carrier_order = demand.loc[
    :, demand.groupby(level=1).sum().sort_values().index
].index.get_level_values(0)

df = demand.unstack().div(1e9)
df = df[df.sum().sort_values().index].loc[carrier_order]

color = n.carriers.set_index("nice_name").color.to_dict()


fig, ax = plt.subplots(figsize=(5, 4))
df.T.plot.bar(stacked=True, color=color, xlabel="", ylabel="Energy demand [PWh]", ax=ax)
ax.legend(title="", frameon=False, bbox_to_anchor=(0, 1.2), loc="upper left")

sns.despine()


fig.savefig(
    "/home/fabian/papers/co2-network/results/baseline/figures/90_nodes/total-demand-bar.png",
    bbox_inches="tight",
)
