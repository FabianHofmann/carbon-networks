#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:14:48 2022

@author: fabian
"""

# Sources and sinks
# Focussing on CO2 network might be enough (turn on/off)
# How to handle import of H2?


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from common import import_network, pypsa_eur_sec, pypsa_eur
from constants import PLOT_SPECS, KINDS

from plotting import (
    create_carrier_network,
    plot_map,
    add_legend,
    get_carrier_network_plotting_data,
)
import geopandas as gpd


sns.set_context("paper")
kind = "carbon"
column = "Optimal Capacity"
alpha = 1
region_alpha = 0.8


n = import_network(
    pypsa_eur_sec
    / "results/your-run-name/postnetworks/elec_s_10_lv1.2__Co2L0-T-H-B-I-A-solar+p3-dist1_2050.nc"
)
regions = gpd.read_file(
    pypsa_eur / "resources/regions_onshore_elec_s_10.geojson"
).set_index("name")

fig, axes = plt.subplots(
    2, 2, figsize=(15, 13), subplot_kw={"projection": ccrs.EqualEarth()}
)

for kind, ax in zip(KINDS, axes.flatten()):

    o = create_carrier_network(n, kind=kind)
    data = get_carrier_network_plotting_data(o, "Optimal Capacity")
    specs = PLOT_SPECS[kind]

    plot_map(
        ax,
        o,
        regions,
        column=column,
        bus_scale=specs["bus_scale"],
        branch_scale=specs["branch_scale"],
        alpha=alpha,
        region_alpha=region_alpha,
        region_cmap=specs["region_cmap"],
    )
    add_legend(
        ax,
        o,
        bus_scale=specs["bus_scale"],
        branch_scale=specs["branch_scale"],
        bus_sizes=specs["bus_sizes"],
        branch_sizes=specs["branch_sizes"],
        alpha=alpha,
    )
    ax.set_extent(regions.total_bounds[[0, 2, 1, 3]])
    ax.set_title(kind.title())

fig.tight_layout()
fig.savefig("map_networks.png", bbox_inches="tight", dpi=300, facecolor="white")
