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
from pypsa.plot import add_legend_circles, add_legend_patches
from matplotlib.gridspec import GridSpec
from common import (
    import_network,
    mock_snakemake,
    get_transmission_links,
    get_carrier_transport,
)
import geopandas as gpd


warnings.filterwarnings("ignore", category=UserWarning)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_operation_map_single",
        run="co2-only",
        clusters=90,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])

n = import_network(snakemake.input.network)
regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
which = "operation"
config = snakemake.config

alpha = 1
region_alpha = 0.6

is_transport = get_transmission_links(n)
transport_carriers = [
    *n.links.carrier[is_transport].unique(),
    *n.lines.carrier.unique(),
]
transport_carriers = n.carriers.nice_name[transport_carriers]

# TODO: ensure this weird adjustment is not necessary anymore
dac = n.links.index[n.links.carrier == "DAC"]
n.links.loc[dac, ["bus0", "bus2"]] = n.links.loc[dac, ["bus2", "bus0"]].values
n.links.loc[dac, ["", "bus2"]] = n.links.loc[dac, ["bus2", "bus0"]].values
n.links_t.p0[dac], n.links_t.p2[dac] = n.links_t.p2[dac], n.links_t.p0[dac].values

s = n.statistics

run = snakemake.wildcards.run
draw_legend = snakemake.params.settings["legend"].get(run, True)

for kind, output in snakemake.output.items():

    # if kind != "hydrogen":
    # continue

    figsize = snakemake.params.settings["figsize"]
    figsize = figsize["single"] if draw_legend else figsize["single_wo_legend"]
    fig, ax = plt.subplots(
        figsize=figsize, subplot_kw={"projection": ccrs.EqualEarth()}
    )

    carriers = config["constants"]["carrier_to_buses"].get(kind, [kind])

    grouper = s.groupers.get_bus_and_carrier
    df = s.dispatch(bus_carrier=carriers, groupby=grouper)
    df = df.drop(transport_carriers, level=2, errors="ignore")
    df = df.rename(lambda x: x.replace(" CC", ""), level=2)
    df = df.groupby(level=[1, 2]).sum().rename(n.buses.location, level=0)

    specs = config["plotting"]["operation_map"][kind]
    bus_scale = float(specs["bus_scale"])
    branch_scale = float(specs["branch_scale"])
    flow_scale = float(specs["flow_scale"])
    legend_kwargs = {"loc": "upper left", "frameon": False}
    unit = specs["unit"]

    bus_sizes = df.sort_index()
    branch_widths = get_carrier_transport(n, kind, config, which)
    flow = pd.concat(branch_widths)
    branch_carriers = get_carrier_transport(n, kind, config, "carrier")
    branch_colors = {
        c: (n.carriers.color.get(vals[0], "lightgrey") if len(vals) else "lightgrey")
        for c, vals in branch_carriers.items()
    }

    n.plot(
        bus_sizes=bus_sizes * bus_scale,
        bus_colors=n.carriers.set_index("nice_name").color,
        bus_alpha=alpha,
        bus_split_circles=True,
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

    if "vmin" in specs:
        vmin = float(specs["vmin"])
    else:
        vmin = regions.price.min()
    if "vmax" in specs:
        vmax = float(specs["vmax"])
    else:
        vmax = regions.price.max()

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

    if draw_legend:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        title = kind.title() if kind != "carbon" else f"Capturing {kind.title()}"
        cbr = fig.colorbar(
            sm,
            ax=ax,
            label=f"Average Price of {title} [{region_unit}]",
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
        ax.set_title(config["labels"].get(kind, kind.title()) + " Operation")

    legend_bus_sizes = specs["bus_sizes"]
    if legend_bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_scale * 1e6 for s in legend_bus_sizes],
            [f"{s} {unit}" for s in legend_bus_sizes],
            legend_kw={"bbox_to_anchor": (0, 1), **legend_kwargs},
        )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )
