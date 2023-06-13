#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:14:48 2022

@author: fabian
"""
import os
import warnings
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pypsa.plot import add_legend_circles, add_legend_patches, add_legend_lines
from common import (
    import_network,
    mock_snakemake,
    get_carrier_consumption,
    get_carrier_transport,
    get_carrier_production,
)

import geopandas as gpd

sns.set_theme(**snakemake.params["theme"])
warnings.filterwarnings("ignore", category=UserWarning)
alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_operation_map",
        kind="carbon",
        design="co2network",
        sequestration=1000,
        ext="png",
    )


n = import_network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
kinds = snakemake.config["constants"]["kinds"]
config = snakemake.config
labels = config["labels"]
which = "operation"
kind = snakemake.wildcards.kind

fig, axes = plt.subplots(
    1,
    2,
    figsize=snakemake.params.settings["figsize"],
    squeeze=False,
    subplot_kw={"projection": ccrs.EqualEarth()},
)

tags = ["production", "consumption"]

for tag, ax in zip(tags, axes.flatten()):
    specs = config["plotting"]["operation_map"][kind]
    bus_scale = float(specs["bus_scale"])
    branch_scale = float(specs["branch_scale"])
    flow_scale = float(specs["flow_scale"])
    legend_kwargs = {"loc": "upper left", "frameon": False}
    unit = specs["region_unit"]

    if tag == "production":
        bus_sizes = get_carrier_production(n, kind, config, which)
    else:
        bus_sizes = get_carrier_consumption(n, kind, config, which)
    branch_widths = get_carrier_transport(n, kind, config, which)
    flow = pd.concat(branch_widths)
    branch_carriers = get_carrier_transport(n, kind, config, "carrier")
    branch_colors = {
        c: (n.carriers.color.get(vals[0], "lightgrey") if len(vals) else "lightgrey")
        for c, vals in branch_carriers.items()
    }

    n.plot(
        bus_sizes=bus_sizes * bus_scale,
        bus_alpha=alpha,
        line_widths=branch_widths["Line"].abs().reindex(n.lines.index, fill_value=0)
        * branch_scale,
        link_widths=branch_widths["Link"].abs().reindex(n.links.index, fill_value=0)
        * branch_scale,
        line_colors=branch_colors["Line"],
        link_colors=branch_colors["Link"],
        flow=flow * flow_scale,
        ax=ax,
        margin=0.2,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
    )

    gen_carriers = n.carriers.loc[bus_sizes.index.unique(1)].sort_values("color")
    add_legend_patches(
        ax,
        gen_carriers.color,
        gen_carriers.nice_name,
        patch_kw={"alpha": alpha},
        legend_kw={"bbox_to_anchor": (0, 0), "ncol": 2, **legend_kwargs},
    )

    ax.set_extent(snakemake.config["plotting"]["extent"])
    if snakemake.params.settings.get("title", True):
        ax.set_title(labels[kind] + " " + tag.title())

    legend_bus_sizes = specs["bus_sizes"]
    if legend_bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_scale * 1e6 for s in legend_bus_sizes],
            [f"{s} {unit}" for s in legend_bus_sizes],
            legend_kw={"bbox_to_anchor": (0, 1), **legend_kwargs},
        )


fig.tight_layout()
fig.savefig(
    snakemake.output.map,
    dpi=300,
)
