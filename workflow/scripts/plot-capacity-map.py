#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:14:48 2022

@author: fabian
"""
import os
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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


column = "Optimal Capacity"
alpha = 1
region_alpha = 0.8
aligned_carriers = ["carbon", "hydrogen"]


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_capacity_map",
        run="full",
        clusters=90,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])

n = import_network(snakemake.input.network)
offshore_regions = gpd.read_file(snakemake.input.offshore_regions).set_index("name")
onshore_regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
config = snakemake.config
which = "capacity"

for kind, output in snakemake.output.items():

    fig, ax = plt.subplots(
        figsize=snakemake.params["settings"]["figsize"],
        subplot_kw={"projection": ccrs.EqualEarth()},
        layout="constrained" if kind in aligned_carriers else None,
    )

    specs = config["plotting"]["capacity_map"][kind]
    unit = specs["unit"]
    unit_conversion = float(specs["unit_conversion"])
    region_unit = specs["region_unit"]
    region_unit_conversion = float(specs["region_unit_conversion"])
    region_cmap = specs["region_cmap"]

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
        bus_sizes=bus_sizes * unit_conversion,
        bus_alpha=alpha,
        line_widths=branch_widths["Line"].reindex(n.lines.index, fill_value=0)
        * branch_scale
        * unit_conversion,
        link_widths=branch_widths["Link"].reindex(n.links.index, fill_value=0)
        * branch_scale
        * unit_conversion,
        line_colors=branch_colors["Line"],
        link_colors=branch_colors["Link"],
        ax=ax,
        margin=0.2,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
    )

    region_data = get_carrier_storage(n, kind, config, which)
    if kind == "carbon":
        offshore_data = region_data.loc[:, "co2 sequestered"]
        onshore_data = region_data.loc[:, "co2 stored"]
        regions = [
            onshore_regions.assign(color=onshore_data),
            offshore_regions.assign(color=offshore_data),
        ]
        regions = pd.concat(regions)
    else:
        region_data = region_data.groupby(level=0).sum()
        regions = onshore_regions.assign(color=region_data)

    regions.color = regions.color * region_unit_conversion

    region_cmap = plt.cm.get_cmap(region_cmap)
    # Create an array of colors from the "Reds" colormap
    cmap = region_cmap(np.linspace(0, 1, 256))
    # Create a new colormap that starts with white and transitions to the colors in "Reds"
    region_cmap = mcolors.LinearSegmentedColormap.from_list(
        "WhiteReds", np.vstack(([1, 1, 1, 1], cmap))
    )

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
            [f"{s} {unit}" for s in legend_bus_sizes],
            legend_kw={"bbox_to_anchor": (0, 1), **legend_kwargs},
        )

    legend_branch_sizes = specs["branch_sizes"]
    if legend_branch_sizes is not None:
        add_legend_lines(
            ax,
            [s * branch_scale for s in legend_branch_sizes],
            [f"{s} {unit}" for s in legend_branch_sizes],
            legend_kw={"bbox_to_anchor": (0, 0.85), "framealpha": 0.7, **legend_kwargs},
        )

    gen_carriers = n.carriers.loc[bus_sizes.index.unique(1)]

    # align figsize for carbon and hydrogen plots
    if kind in aligned_carriers:
        legend_kw_overwrite = {"bbox_to_anchor": (0, 1.25), "loc": "upper left"}
        if kind == "hydrogen":
            legend_kw_overwrite["ncol"] = 1
    else:
        legend_kw_overwrite = {}

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
            **legend_kw_overwrite,
        },
    )

    # fig.set_title(kind.title() + " Production and Storage Capacity Map")

    fig.savefig(
        output,
        # bbox_inches="tight",
        dpi=300,
        facecolor="white",
    )
