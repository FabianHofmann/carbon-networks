import pandas as pd
import numpy as np
import pypsa
from pypsa.descriptors import Dict
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches
from common import get_transmission_links, get_generation_carriers
from cartopy import crs as ccrs


def plot_map(
    ax,
    n,
    regions,
    data,
    bus_scale=1.5e-5,
    branch_scale=3e-3,
    alpha=1,
    branch_alpha=1,
    region_alpha=0.8,
    region_cmap="Blues",
    region_unit="GWh",
    **kwargs,
):

    n.bus_scale = bus_scale
    n.branch_scale = branch_scale
    n.alpha = alpha

    # with plt.rc_context({"patch.linewidth": 0.0}):
    n.plot(
        bus_sizes=data.bus_sizes * bus_scale,
        bus_alpha=alpha,
        line_widths=data.line_widths * branch_scale,
        link_widths=data.link_widths * branch_scale,
        link_colors=data.branch_color,
        line_colors=data.branch_color,
        line_alpha=branch_alpha,
        link_alpha=branch_alpha,
        flow=data.flow * branch_scale if data.flow is not None else None,
        ax=ax,
        margin=0.2,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
        **kwargs,
    )
    if data.region_color is not None:
        regions = regions.assign(color=data.region_color / 1e6)
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
                "label": n.kind.title() + f" storage capacity [{region_unit}]",
                "orientation": "horizontal",
                "shrink": 0.8,
                "pad": 0.05,
                "aspect": 30,
                "alpha": region_alpha,
            },
        )


def add_legend(
    ax,
    n,
    bus_sizes=[10_000, 5_000],
    branch_sizes=[2_000, 1_000],
    bus_scale=None,
    branch_scale=None,
    alpha=None,
    gen_carriers=None,
    **kwargs,
):
    legend_kwargs = {"loc": "upper left", "frameon": False, **kwargs}

    if bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_scale for s in bus_sizes],
            [f"{s // 1000} GW" for s in bus_sizes],
            legend_kw={"bbox_to_anchor": (1, 1), **legend_kwargs},
        )
    if branch_sizes is not None:
        add_legend_lines(
            ax,
            [s * branch_scale for s in branch_sizes],
            [f"{s // 1000} GW" for s in branch_sizes],
            legend_kw={"bbox_to_anchor": (1, 0.9), **legend_kwargs},
        )

    if gen_carriers is None:
        gen_carriers = n.carriers.loc[get_generation_carriers(n)]
    add_legend_patches(
        ax,
        gen_carriers.color,
        gen_carriers.nice_name,
        patch_kw={"alpha": alpha},
        legend_kw={"bbox_to_anchor": (1, 0.1), **legend_kwargs, "loc": "lower left"},
    )
    branch_carriers = n.carriers.loc[
        n.links[get_transmission_links(n)].carrier.unique()
    ]
    add_legend_patches(
        ax,
        branch_carriers.color,
        branch_carriers.nice_name,
        patch_kw={"alpha": alpha},
        legend_kw={"bbox_to_anchor": (1, 0.0), **legend_kwargs, "loc": "lower left"},
    )
