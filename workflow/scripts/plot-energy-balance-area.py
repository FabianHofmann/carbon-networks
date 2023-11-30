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
import matplotlib.dates as mdates

from common import (
    sort_rows_by_relative_diff,
    import_network,
    mock_snakemake,
)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_energy_balance_area",
        ext="pdf",
        clusters=90,
        run="co2-only",
    )

sns.set_theme(**snakemake.params["theme"])
plt.rc("patch", linewidth=0)

config = snakemake.config
kind_to_carrier = {"hydrogen": "H2", "carbon": "co2 stored", "electricity": "AC"}


n = import_network(snakemake.input.network)
key = snakemake.params.label
balance = n.statistics.energy_balance(aggregate_time=False)


for kind, output in snakemake.output.items():

    carrier = kind_to_carrier.get(kind, kind)
    label = config["labels"].get(kind, kind.title())

    norm = 1e6
    unit = "TWh" if kind not in ["carbon", "co2"] else "Mt"

    ds = balance.xs(carrier, level=2)
    ds = ds.div(norm)
    ds = ds.droplevel(0).dropna(how="all").fillna(0).groupby(level=0).sum()
    # ds = ds.sort_values(ds.columns[0], ascending=False)
    ds = sort_rows_by_relative_diff(ds)

    fig, ax = plt.subplots(1, 1, figsize=snakemake.params.settings["figsize"])

    colors = n.carriers.set_index("nice_name").color.to_dict()
    ylim = ds.where(ds >= 0).sum().max() * 1.1
    pos = ds[ds.gt(0).any(axis=1)]
    neg = ds[ds.le(0).any(axis=1)]

    pos[::-1].where(lambda ds: ds > 0).T.plot.area(
        color=colors, ax=ax, stacked=True, alpha=0.8, lw=0, rot=0, ylim=(-ylim, ylim)
    )
    neg.where(lambda ds: ds <= 0).T.plot.area(
        color=colors, ax=ax, stacked=True, alpha=0.8, lw=0, rot=0, ylim=(-ylim, ylim)
    )

    ax.axhline(0, color="k", lw=1)
    ax.set_ylabel(f"{label} [{unit}]")
    ax.set_xlim(ds.columns[0], ds.columns[-1])
    ax.grid(axis="y", alpha=0.5)
    if snakemake.params.settings.get("title", True):
        ax.set_title(f"{label} Balance")
    ax.set_xlabel("")

    formatter = mdates.DateFormatter("%b")  # %b gives the abbreviated month name
    ax.xaxis.set_major_formatter(formatter)
    fig.autofmt_xdate(ha="left", rotation=0)

    h, l = ax.get_legend_handles_labels()
    order = pd.Index([*pos.index, *neg.index]).drop_duplicates(keep="first")
    legend = pd.Series(h, index=l)[lambda ds: ~ds.index.duplicated()][order]
    ax.legend(
        legend.values,
        legend.index,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        frameon=False,
        ncol=1,
    )

    sns.despine()
    fig.savefig(output, dpi=300, bbox_inches="tight")
