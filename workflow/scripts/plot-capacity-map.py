#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:14:48 2022

@author: fabian
"""
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pypsa.plot import add_legend_circles, add_legend_patches, add_legend_lines
from common import (
    import_network,
    mock_snakemake,
    get_carrier_storage,
    get_carrier_transport,
    get_carrier_production,
)

import geopandas as gpd

sns.set_theme(style="white", context="paper", rc={"patch.linewidth": 0.1}, font="serif")
column = "Optimal Capacity"
alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_capacity_map",
        kind="carbon",
        design="co2network",
        sequestration=200,
    )


n = import_network(snakemake.input.network)
offshore_regions = gpd.read_file(snakemake.input.offshore_regions).set_index("name")
onshore_regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
kind = snakemake.wildcards.kind
config = snakemake.config
which = "capacity"


fig, ax = plt.subplots(
    figsize=(7, 6),
    subplot_kw={"projection": ccrs.EqualEarth()},
)


regions = offshore_regions if kind == "carbon" else onshore_regions
unit = "GW" if kind not in ["carbon", "co2"] else "kt/h"

specs = config["plotting"]["capacity_map"][kind]
bus_scale = float(specs["bus_scale"])
branch_scale = float(specs["branch_scale"])

bus_sizes = get_carrier_production(n, kind, config, which) * bus_scale
branch_widths = get_carrier_transport(n, kind, config, which)
branch_carriers = get_carrier_transport(n, kind, config, "carrier")
branch_colors = {
    c: (n.carriers.color.get(vals[0], "lightgrey") if len(vals) else "lightgrey")
    for c, vals in branch_carriers.items()
}

n.plot(
    bus_sizes=bus_sizes,
    bus_alpha=alpha,
    line_widths=branch_widths["Line"].reindex(n.lines.index, fill_value=0)
    * branch_scale,
    link_widths=branch_widths["Link"].reindex(n.links.index, fill_value=0)
    * branch_scale,
    line_colors=branch_colors["Line"],
    link_colors=branch_colors["Link"],
    ax=ax,
    margin=0.2,
    color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
    geomap="10m",
)

region_data = get_carrier_storage(n, kind, config, which)
region_data = region_data.groupby(level=0).sum()
region_cmap = specs["region_cmap"]
region_unit = specs["region_unit"]

regions = regions.assign(color=region_data / 1e6)
regions.plot(
    ax=ax,
    column="color",
    cmap=region_cmap,
    vmin=0,
    alpha=region_alpha,
    facecolor="whitesmoke",
    edgecolor="black",
    aspect="equal",
    transform=ccrs.PlateCarree(),
    linewidth=0.0,
    legend=True,
    legend_kwds={
        "label": kind.title() + f" storage capacity [{region_unit}]",
        "orientation": "horizontal",
        "shrink": 0.8,
        "pad": 0.05,
        "aspect": 30,
        "alpha": region_alpha,
    },
)
ax.set_extent(snakemake.config["plotting"]["extent"])

legend_bus_sizes = specs["bus_sizes"]
legend_kwargs = {"framealpha": 0.7, "frameon": True, "loc": "upper left"}
if legend_bus_sizes is not None:
    add_legend_circles(
        ax,
        [s * bus_scale for s in legend_bus_sizes],
        [f"{s // 1000} {unit}" for s in legend_bus_sizes],
        legend_kw={"bbox_to_anchor": (0, 1), **legend_kwargs},
    )

legend_branch_sizes = specs["branch_sizes"]
if legend_branch_sizes is not None:
    add_legend_lines(
        ax,
        [s * branch_scale for s in legend_branch_sizes],
        [f"{s // 1000} {unit}" for s in legend_branch_sizes],
        legend_kw={"bbox_to_anchor": (0, 0.85), "framealpha": 0.7, **legend_kwargs},
    )

gen_carriers = n.carriers.loc[bus_sizes.index.unique(1)]
add_legend_patches(
    ax,
    gen_carriers.color,
    gen_carriers.nice_name,
    patch_kw={"alpha": alpha},
    legend_kw={
        "bbox_to_anchor": (0, 1),
        "ncol": 2,
        "frameon": False,
        "loc": "lower left",
    },
)


# fig.set_title(kind.title() + " Production and Storage Capacity Map")

fig.savefig(
    snakemake.output.map,
    bbox_inches="tight",
    dpi=300,
    facecolor="white",
)
