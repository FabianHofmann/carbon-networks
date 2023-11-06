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
        "plot_sankey_diagramm", ext="png", clusters=90, run="co2-only"
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
    "gas for industry": 0.2,
    "methanol": 0.248009559,
    "process emissions": 1,
    # "solid biomass": 0, # is net neutral
    # "biogas": 0.0, # is net neutral
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
drop = drop[~drop.get_level_values(2).str.contains("pipeline")]
# test for net neutrality
# assert source.loc[drop].round(0).equals(target.loc[drop].round(0))
source = source.drop(drop)
target = target.drop(drop)

# test for net neutrality
# assert round(target.sum(), 0) == round(source.sum(), 0)


target = target.reset_index(0, name="value").rename(columns={"level_0": "target"})
source = source.reset_index(0, name="value").rename(columns={"level_0": "source"})

data = target.assign(source=source.source.reindex(target.index)).dropna(
    subset=["source"]
)
data = data.droplevel(0)
data = data[data.value > 10]


int_map = {v: k for k, v in enumerate(set([*data.source, *data.target]))}
labels = n.carriers.nice_name[int_map.keys()].str.replace("$_2$", "₂")
colors = n.carriers.color[int_map.keys()]
link_color = n.carriers.color[data.value.index]
link_label = n.carriers.nice_name[data.index].str.replace("$_2$", "₂")
link_label["oil emissions"] = "Aviation & Petrochemicals"
link_meta = pd.concat([link_color, link_label], axis=1)
link_meta.drop_duplicates(inplace=True)
link_meta.sort_values(by="nice_name", inplace=True)

# Define the links
source = data.source.replace(int_map)
target = data.target.replace(int_map)
value = data.value

# Create the Sankey diagram
# Create a new list of labels with a white color
shadow_labels = [
    '<span style="color:black; text-shadow: -1px 0 white, 0 1px white, 1px 0 white, 0 -1px white;">{}</span>'.format(
        label
    )
    for label in labels
]


fig = go.Figure(
    data=go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=shadow_labels,
            color=colors,
        ),
        link=dict(source=source, target=target, value=value, color=link_color),
        legend="legend",
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

# Update the layout
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
    ),
    margin=dict(l=20, r=20, t=20, b=20),
)

fig.write_image(snakemake.output[0], width=800, height=600)
