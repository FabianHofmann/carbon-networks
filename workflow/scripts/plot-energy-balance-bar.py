#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 18:11:52 2023

@author: fabian
"""
import os
import pandas as pd
import pypsa
import textwrap
import matplotlib.pyplot as plt
import seaborn as sns
from common import (
    import_network,
    mock_snakemake,
    sort_rows_by_diff,
    get_ordered_handles_labels,
)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_energy_balance_bar",
        ext="png",
        clusters=90,
        comparison="default",
    )

sns.set_theme(**snakemake.params["theme"])

config = snakemake.config
norm = 1e6
wrapper = textwrap.TextWrapper(width=25)


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    balance = n.statistics.energy_balance()

    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    df[key] = balance


df = pd.concat(df, axis=1)

for kind, output in snakemake.output.items():

    kind_to_carrier = {"hydrogen": "H2", "carbon": "co2 stored", "electricity": "AC"}
    carrier = kind_to_carrier.get(kind, kind)
    label = config["labels"].get(kind, kind.title())
    unit = "TWh" if kind not in ["carbon", "co2"] else "Mt"

    ds = df.xs(carrier, level=2)
    ds = ds.div(norm)
    ds = ds[ds.round(2) != 0].droplevel(0).dropna(how="all").fillna(0)
    ds = ds.groupby(level=0).sum()
    ds = ds.sort_values(ds.columns[0], ascending=False)
    if kind == "co2":
        ds.drop("CO$_2$", inplace=True, errors="ignore")
    ds = sort_rows_by_diff(ds)

    fig, ax = plt.subplots(
        1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
    )

    colors = n.carriers.set_index("nice_name").color.to_dict()
    ds.T.plot.bar(color=colors, ax=ax, stacked=True, alpha=0.8, lw=0, rot=90)

    ax.axhline(0, color="k", lw=1)
    ax.set_ylabel(f"{label} [{unit}]")
    ax.grid(axis="y", alpha=0.5)
    if snakemake.params.settings.get("title", True):
        ax.set_title(f"{label} Balance")

    handles, labels = get_ordered_handles_labels(ax, ds, wrap=25)
    ax.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        frameon=False,
        ncol=1,
    )

    sns.despine()
    fig.savefig(output, dpi=300)
