#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:14:48 2022

@author: fabian
"""
import os
import warnings
import textwrap
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
wrapper = textwrap.TextWrapper(width=18)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_balance_map_dedicated",
        clusters=90,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])
plt.rc("patch", linewidth=0.1)

onregions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
offregions = gpd.read_file(snakemake.input.offshore_regions).set_index("name")
networks = [
    import_network(
        path,
        revert_dac=True,
        offshore_sequestration=True,
        offshore_regions=offregions,
    )
    for path in snakemake.input.networks
]
config = snakemake.config
labels = config["labels"]
alpha = config["plotting"]["balance_map"]["alpha"]
region_alpha = config["plotting"]["balance_map"]["region_alpha"]

fig, axes = plt.subplots(
    2,
    2,
    figsize=(8, 8),
    subplot_kw={"projection": ccrs.EqualEarth()},
    layout="constrained",
)

bounds = snakemake.config["plotting"]["extent"]
bounds = [-10, 28, 36, 70]  # adjust for better alignment with latex document

for n, axs, right_subplot in zip(networks, axes.T, [False, True]):
    s = n.statistics
    colors = n.carriers.set_index("nice_name").color
    run = n.meta["wildcards"]["run"]

    for kind, ax in zip(["carbon", "hydrogen"], axs):
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
            "frameon": False,
            "framealpha": 1,
            "edgecolor": "None",
        }
        unit = specs["unit"]

        bus_sizes = df.sort_index()
        flow = s.transmission(groupby=False, bus_carrier=carriers)
        branch_colors = {c: colors[carrier] for c, carrier in transmission_carriers}
        fallback = pd.Series()

        # plot sequestration sinks as full circles, watch out in current pypsa verion
        # the bus sizes are reduced by factor 2 if split circles is activated!
        # https://github.com/PyPSA/PyPSA/issues/799
        bus_sizes = bus_sizes * 2

        if kind == "carbon":
            sequestration_sizes = -bus_sizes.loc[:, ["CO$_2$ Sequestration"]] / 2
            bus_sizes = bus_sizes.drop("CO$_2$ Sequestration", level=1)
            n.plot(
                bus_sizes=sequestration_sizes * bus_scale,
                bus_colors=colors,
                bus_alpha=alpha,
                line_widths=0,
                link_widths=0,
                ax=ax,
                color_geomap=False,
                geomap=True,
                boundaries=bounds,
            )
            offregions.plot(
                ax=ax,
                facecolor="None",
                edgecolor="darkgrey",
                linewidth=0.1,
                alpha=region_alpha,
                transform=ccrs.PlateCarree(),
                aspect="equal",
            )

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
            boundaries=bounds,
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

        onregions["price"] = price.reindex(onregions.index).fillna(0)
        region_unit = specs["region_unit"]
        cmap = specs["region_cmap"]
        vmin, vmax = specs["vmin"], specs["vmax"]

        onregions.plot(
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

        if right_subplot:
            sm = plt.cm.ScalarMappable(
                cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax)
            )
            title = kind.title() if kind != "carbon" else f"Capturing {kind.title()}"
            cbr = fig.colorbar(
                sm,
                ax=axes[0] if kind == "carbon" else axes[1],
                label=f"Average Marginal Price of {title} [{region_unit}]",
                shrink=0.6,
                pad=0.05,
                aspect=50,
                alpha=region_alpha,
                orientation="horizontal",
            )
            cbr.outline.set_edgecolor("None")

            pad = 0.18
            prod_carriers = (
                bus_sizes[bus_sizes > 0].index.unique("carrier").sort_values()
            )
            cons_carriers = (
                bus_sizes[bus_sizes < 0]
                .index.unique("carrier")
                .difference(prod_carriers)
                .sort_values()
            )

            # fix bug related to falsely clipped normalization
            if kind == "hydrogen" and "H$_2$ For Industry" in prod_carriers:
                prod_carriers = prod_carriers.difference(["H$_2$ For Industry"])
                cons_carriers = cons_carriers.union(["H$_2$ For Industry"])
            elif kind == "carbon" and "CO$_2$ Sequestration" not in cons_carriers:
                cons_carriers = cons_carriers.union(["CO$_2$ Sequestration"])

            add_legend_patches(
                ax,
                colors[prod_carriers],
                prod_carriers.map(wrapper.fill),
                patch_kw={"alpha": alpha},
                legend_kw={
                    "loc": "upper left",
                    "bbox_to_anchor": (1, 1.1),
                    "title": "Production",
                    **legend_kwargs,
                },
            )

            add_legend_patches(
                ax,
                colors[cons_carriers],
                cons_carriers.map(wrapper.fill),
                patch_kw={"alpha": alpha},
                legend_kw={
                    "loc": "upper left",
                    "bbox_to_anchor": (1, 0.56),
                    "title": "Consumption",
                    **legend_kwargs,
                },
            )

            # only use one legend for both branches and flows
            legend_bus_sizes = specs["bus_sizes"][:1]
            add_legend_circles(
                ax,
                [s * bus_scale * 1e6 for s in legend_bus_sizes],
                [f"{s} {unit}" for s in legend_bus_sizes],
                legend_kw={
                    "loc": "lower left",
                    "bbox_to_anchor": (1, 0.07),
                    **legend_kwargs,
                },
            )
            # only use one legend for both branches and flows
            legend_branch_sizes = specs["branch_sizes"][:1]
            if legend_branch_sizes is not None:
                add_legend_lines(
                    ax,
                    [s * branch_scale * 1e6 for s in legend_branch_sizes],
                    [f"{s} {unit}" for s in legend_branch_sizes],
                    legend_kw={
                        "loc": "lower left",
                        "bbox_to_anchor": (1, 0),
                        **legend_kwargs,
                    },
                )

        ax.set_extent(bounds)
        if snakemake.params.settings.get("title", True):
            carrier_label = labels.get(kind, kind.title())
            title = f"{carrier_label} Balance ({labels[run]} Model)"
            ax.set_title(title)

fig.savefig(
    snakemake.output.figure,
    dpi=350,
)
