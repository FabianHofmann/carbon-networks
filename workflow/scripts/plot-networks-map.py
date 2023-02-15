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
from common import import_network, mock_snakemake

from plotting import (
    create_carrier_network,
    plot_map,
    add_legend,
    get_carrier_network_plotting_data,
)
import geopandas as gpd


sns.set_context("paper")
column = "Optimal Capacity"
alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_networks_map",
        simpl="",
        lv=1.2,
        clusters=181,
        opts="",
        sector_opts="Co2L0-365H-T-H-B-I-A-solar+p3-linemaxext15-seq200",
        planning_horizons=2050,
    )


n = import_network(snakemake.input.network)
offshore_regions = gpd.read_file(snakemake.input.offshore_regions).set_index("name")
onshore_regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
kinds = snakemake.config["constants"]["kinds"]

fig, axes = plt.subplots(
    2,
    2,
    figsize=(15, 13),
    subplot_kw={"projection": ccrs.EqualEarth()},
)

for kind, ax in zip(kinds, axes.flatten()):

    carriers = snakemake.config["constants"]["carrier_to_buses"][kind]
    o = create_carrier_network(n, kind, carriers)
    data = get_carrier_network_plotting_data(o, "Optimal Capacity")
    specs = snakemake.config["plotting"]["capacity_map"][kind]

    if kind == "carbon":
        regions = offshore_regions
    else:
        regions = onshore_regions

    plot_map(
        ax,
        o,
        regions,
        data,
        bus_scale=float(specs["bus_scale"]),
        branch_scale=float(specs["branch_scale"]),
        alpha=alpha,
        branch_alpha=alpha,
        region_alpha=region_alpha,
        region_cmap=specs["region_cmap"],
        region_unit=specs["region_unit"],
    )
    fig.canvas.draw()
    add_legend(
        ax,
        o,
        bus_scale=float(specs["bus_scale"]),
        branch_scale=float(specs["branch_scale"]),
        bus_sizes=specs["bus_sizes"],
        branch_sizes=specs["branch_sizes"],
        alpha=alpha,
    )
    ax.set_extent(regions.total_bounds[[0, 2, 1, 3]])
    ax.set_title(kind.title())

fig.tight_layout()
fig.savefig(
    snakemake.output.map,
    bbox_inches="tight",
    dpi=300,
    facecolor="white",
)
