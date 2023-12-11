#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 10:55:25 2023

@author: fabian
"""


import plotly.graph_objects as go
import plotly.io as pio
import os
import pandas as pd
from common import import_network, mock_snakemake


pio.renderers.default = "svg"

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_sankey_diagramm", ext="png", clusters=90, run="net-negative-0.1/full"
    )


n = import_network(snakemake.input.network)
n.carriers.loc["co2", "nice_name"] = "CO$_2$ Atmospheric"
s = n.statistics


carbon_content = {
    "co2 stored": 1,
    "gas": 0.2,
    "co2": 1,
    "co2 sequestered": 1,
    "oil": 0.26,
    "gas for industry": 0,
    "methanol": 0.248009559,
    "process emissions": 1,
    # "solid biomass": 0, # is net neutral
    # "biogas": 0., # is net neutral
    # "solid biomass for industry": 0, # is net neutral
}

supply, withdrawal = [], []
for k, v in carbon_content.items():
    supply.append(s.supply(bus_carrier=k, nice_names=False) * v)
    withdrawal.append(s.withdrawal(bus_carrier=k, nice_names=False) * v)
withdrawal = pd.concat(withdrawal, keys=list(carbon_content.keys()))
supply = pd.concat(supply, keys=list(carbon_content.keys()))

rename = {c: c[:-3] for c in n.carriers.index if c.endswith(" CC")}


source = withdrawal.rename(rename, level=2).groupby(level=[0, 1, 2]).sum()
target = supply.rename(rename, level=2).groupby(level=[0, 1, 2]).sum()

rename = {
    "agriculture machinery oil": "agriculture machinery oil emissions",
    "kerosene for aviation": "oil emissions",
    "naphtha for industry": "oil emissions",
    "shipping methanol": "shipping methanol emissions",
}
source = source.rename(index=rename, level=2).groupby(level=[0, 1, 2]).sum()


drop = source.index.intersection(target.index)
# drop = drop[~drop.get_level_values(2).str.contains("pipeline")]
# test for net neutrality
# assert source.loc[drop].round(0).equals(target.loc[drop].round(0))
source = source.drop(drop)
target = target.drop(drop)

# test for net neutrality
# assert round(target.sum(), 0) == round(source.sum(), 0)


target = target.reset_index(0, name="value").rename(columns={"level_0": "target"})
source = source.reset_index(0, name="value").rename(columns={"level_0": "source"})

data = target.assign(source=source.source.reindex(target.index))
data["source"] = data.source.fillna("fossil " + data.target.loc[["Generator"]])
data = data.droplevel(0).dropna(subset=["source"])
data = data[data.value > 10]
data = data.sort_values(by="value", ascending=True)

carriers = n.carriers.copy()
fossils = data.source[data.source.str.startswith("fossil")]
fossils_map = fossils.str.replace("fossil ", "")
carriers = pd.concat([carriers, carriers.loc[fossils_map].rename(index=fossils)])
carriers.loc[fossils, "nice_name"] = "Fossil " + carriers.loc[fossils, "nice_name"]

int_map = {v: k for k, v in enumerate(set([*data.source, *data.target]))}
labels = carriers.nice_name[int_map.keys()].str.replace("$_2$", "₂")
colors = carriers.color[int_map.keys()]
link_color = carriers.color[data.value.index]
link_label = carriers.nice_name[data.index].str.replace("$_2$", "₂")
link_label["oil emissions"] = "Aviation & Petrochemicals"
link_meta = pd.concat([link_color, link_label], axis=1)
link_meta.drop_duplicates(inplace=True)
link_meta.sort_values(by="nice_name", inplace=True)
# Define the links

source = data.source.replace(int_map)
target = data.target.replace(int_map)
value = data.value.round(0).astype(int)

# Create the Sankey diagram
fig = go.Figure(
    data=go.Sankey(
        arrangement="freeform",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color=colors,
        ),
        link=dict(source=source, target=target, value=value, color=link_color),
    )
)

for label, color in zip(link_meta.nice_name, link_meta.color):
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(size=10, color=color),
            showlegend=True,
            name=label,
        )
    )

fig.update_layout(
    font=dict(size=12),
    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    plot_bgcolor="white",
    legend=dict(
        orientation="h",
        yanchor="top",
        y=-0.2,
        xanchor="center",
        x=0.5,
        traceorder="normal",
        itemsizing="constant",
        bgcolor="rgba(255,255,255,0.8)",
    ),
    margin=dict(l=30, r=30, t=30, b=30),
    autosize=False,
    width=800,
    height=500,
)


fig.show()
fig.write_image(snakemake.output[0], scale=3)
