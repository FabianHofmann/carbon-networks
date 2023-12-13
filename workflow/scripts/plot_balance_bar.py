#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 18:11:52 2023

@author: fabian
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import (
    import_network,
    mock_snakemake,
    sort_rows_by_relative_diff,
    get_ordered_handles_labels,
)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_balance_bar",
        ext="png",
        clusters=40,
        comparison="emission-reduction-full",
    )

sns.set_theme(**snakemake.params["theme"])

config = snakemake.config

df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    balance = n.statistics.energy_balance()

    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    df[key] = balance


df = pd.concat(df, axis=1).fillna(0)
df = df.rename(lambda k: k.replace(" CC", ""), level=1).groupby(level=[0, 1, 2]).sum()

for kind, output in snakemake.output.items():
    if kind.startswith("table"):
        continue
    kind_to_carrier = {"hydrogen": "H2", "carbon": "co2 stored", "electricity": "AC"}
    carrier = kind_to_carrier.get(kind, kind)
    label = config["labels"].get(kind, kind.title())

    if kind in ["carbon", "co2"]:
        unit = "Mt/a"
        norm = 1e6
        fmt = ".0f"
    elif kind in ["electricity", "hydrogen", "heat", "oil"]:
        unit = "PWh"
        norm = 1e9
        fmt = ".2f"
    else:
        unit = "TWh"
        norm = 1e6
        fmt = ".0f"

    ds = df.xs(carrier, level=2)
    ds = ds.div(norm).droplevel(0)
    reduced = ds.where(ds.round(2) != 0).dropna(how="all").fillna(0)
    if not reduced.empty:
        ds = reduced
    ds = ds.groupby(level=0).sum()
    ds = ds.sort_values(ds.columns[0], ascending=False)
    if kind == "co2":
        ds.drop("CO$_2$", inplace=True, errors="ignore")
    ds = sort_rows_by_relative_diff(ds)

    fig, ax = plt.subplots(
        1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
    )

    colors = n.carriers.set_index("nice_name").color.to_dict()
    ds.T.plot.bar(color=colors, ax=ax, stacked=True, alpha=0.8, rot=90)
    total = ds[ds > 0].sum()
    pad = ax.get_ylim()[1] * 0.01
    for i, val in enumerate(total):
        ax.text(
            i,
            val + pad,
            f"{val:{fmt}}",
            ha="center",
            va="bottom",
            fontsize=8,
            color="grey",
        )

    ax.axhline(0, color="k", lw=1)
    ax.set_ylabel(f"{label} [{unit}]")
    # ax.grid(axis="y", alpha=0.5)
    if snakemake.params.settings.get("title", True):
        ax.set_title(f"{label} Balance")

    handles, labels = get_ordered_handles_labels(ax, ds, wrap=22)
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
    ds.round(3).to_csv(snakemake.output[f"table-{kind}"])
