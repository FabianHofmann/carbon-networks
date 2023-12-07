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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pypsa.plot import add_legend_circles, add_legend_patches, add_legend_lines
from pypsa.statistics import get_transmission_carriers
from common import (
    import_network,
    mock_snakemake,
)
import geopandas as gpd


warnings.filterwarnings("ignore", category=UserWarning)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_balance_map",
        run="full",
        clusters=90,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])
plt.rc("patch", linewidth=0.1)

n = import_network(snakemake.input.network, revert_dac=True)
regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
config = snakemake.config
labels = config["labels"]

alpha = config["plotting"]["balance_map"]["alpha"]
region_alpha = config["plotting"]["balance_map"]["region_alpha"]

s = n.statistics
colors = n.carriers.set_index("nice_name").color

run = snakemake.wildcards.run

for kind, output in snakemake.output.items():

    fig, ax = plt.subplots(
        figsize=snakemake.params.settings["figsize"],
        subplot_kw={"projection": ccrs.EqualEarth()},
        layout="constrained",
    )

    carriers = config["constants"]["carrier_to_buses"].get(kind, [kind])

    grouper = s.groupers.get_bus_and_carrier
    df = s.dispatch(bus_carrier=carriers, groupby=grouper)
    _ = get_transmission_carriers(n, bus_carrier=carriers)
    transmission_carriers = _.set_levels(
        n.carriers.nice_name[_.get_level_values(1)], level=1
    )
    sub = df.loc[["Link"]].drop(
        transmission_carriers.unique(1), level=2, errors="ignore"
    )
    df = pd.concat([df.drop("Link"), sub])
    df = df.rename(lambda x: x.replace(" CC", ""), level=2)
    df = df.groupby(level=[1, 2]).sum().rename(n.buses.location, level=0)
    df = df[df.abs() > 1]

    specs = config["plotting"]["balance_map"][kind]
    bus_scale = float(specs["bus_scale"])
    branch_scale = float(specs["branch_scale"])
    flow_scale = float(specs["flow_scale"])
    legend_kwargs = {
        "loc": "upper left",
        "frameon": True,
        "framealpha": 1,
        "edgecolor": "None",
    }
    unit = specs["unit"]

    bus_sizes = df.sort_index()
    flow = s.transmission(groupby=False, bus_carrier=carriers)
    branch_colors = {c: colors[carrier] for c, carrier in transmission_carriers}
    fallback = pd.Series()

    n.plot(
        bus_sizes=bus_sizes * bus_scale,
        bus_colors=colors,
        bus_alpha=alpha,
        bus_split_circles=True,
        line_widths=flow.get("Line", fallback)
        .abs()
        .reindex(n.lines.index, fill_value=0)
        * branch_scale,
        link_widths=flow.get("Link", fallback)
        .abs()
        .reindex(n.links.index, fill_value=0)
        * branch_scale,
        line_colors=branch_colors.get("Line", "lightgrey"),
        link_colors=branch_colors.get("Link", "lightgrey"),
        flow=flow * flow_scale,
        ax=ax,
        margin=0.2,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
        boundaries=snakemake.config["plotting"]["extent"],
    )

    buses = n.buses.query("carrier in @carriers").index
    price = (
        n.buses_t.marginal_price.mean()
        .reindex(buses)
        .rename(n.buses.location)
        .groupby(level=0)
        .mean()
    )
    if kind == "carbon":
        price = price - n.global_constraints.mu["CO2Limit"]

    regions["price"] = price.reindex(regions.index).fillna(0)
    region_unit = specs["region_unit"]
    cmap = specs["region_cmap"]

    if set(["vmin", "vmax"]).issubset(specs):
        vmin, vmax = specs["vmin"], specs["vmax"]
    else:
        vmin, vmax = regions.price.min(), regions.price.max()

    regions.plot(
        ax=ax,
        column="price",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        edgecolor="None",
        linewidth=0,
        alpha=region_alpha,
        transform=ccrs.PlateCarree(),
        aspect="equal",
    )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    title = kind.title() if kind != "carbon" else f"Capturing {kind.title()}"
    cbr = fig.colorbar(
        sm,
        ax=ax,
        label=f"Average Marginal Price of {title} [{region_unit}]",
        shrink=1,
        pad=0.03,
        aspect=50,
        alpha=region_alpha,
        orientation="horizontal",
    )
    cbr.outline.set_edgecolor("None")

    pad = 0.18
    carriers = n.carriers.set_index("nice_name")
    prod_carriers = bus_sizes[bus_sizes > 0].index.unique("carrier").sort_values()
    cons_carriers = (
        bus_sizes[bus_sizes < 0]
        .index.unique("carrier")
        .difference(prod_carriers)
        .sort_values()
    )

    add_legend_patches(
        ax,
        carriers.color[prod_carriers],
        prod_carriers,
        patch_kw={"alpha": alpha},
        legend_kw={
            "bbox_to_anchor": (0, -pad),
            "ncol": 1,
            "title": "Production",
            **legend_kwargs,
        },
    )

    add_legend_patches(
        ax,
        carriers.color[cons_carriers],
        cons_carriers,
        patch_kw={"alpha": alpha},
        legend_kw={
            "bbox_to_anchor": (0.5, -pad),
            "ncol": 1,
            "title": "Consumption",
            **legend_kwargs,
        },
    )

    ax.set_extent(snakemake.config["plotting"]["extent"])
    if snakemake.params.settings.get("title", True):
        carrier_label = labels.get(kind, kind.title())
        title = f"{carrier_label} Balance ({labels[run]} Model)"
        ax.set_title(title)

    legend_bus_sizes = specs["bus_sizes"]
    if legend_bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_scale * 1e6 for s in legend_bus_sizes],
            [f"{s} {unit}" for s in legend_bus_sizes],
            legend_kw={"bbox_to_anchor": (0, 1), **legend_kwargs},
        )
    legend_branch_sizes = specs["branch_sizes"]
    if legend_branch_sizes is not None:
        add_legend_lines(
            ax,
            [s * branch_scale * 1e6 for s in legend_branch_sizes],
            [f"{s} {unit}" for s in legend_branch_sizes],
            legend_kw={"bbox_to_anchor": (0, 0.85), **legend_kwargs},
        )

    fig.savefig(
        output,
        dpi=300,
    )
