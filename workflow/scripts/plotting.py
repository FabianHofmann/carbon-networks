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


def create_carrier_network(n, kind, carriers, include_eu=False):
    """Create a network with only buses and components of a certain kind.

    Producing technologies are represented by generators and single-linked links.
    Withdrawing technologies are represented by loads.
    """
    buses = n.buses.query("carrier in @carriers").index

    o = pypsa.Network(snapshots=n.snapshots)
    o.madd("Bus", n.buses.index, **n.buses, marginal_price=n.buses_t.marginal_price)
    o.buses.loc[buses, "kind"] = True
    o.buses.kind.fillna(False, inplace=True)

    for c in n.one_port_components:
        if n.df(c).empty:
            continue

        df = n.df(c).query("bus in @buses")
        o.madd(c, df.index, **df, p=n.pnl(c).p[df.index])

    for c in n.branch_components:

        if n.df(c).empty:
            continue

        df = n.df(c)

        if not include_eu:
            df = df[df.bus0.map(n.buses.location).ne("EU")]

        bus_cols = [c for c in df.columns if c.startswith("bus")]
        eff_cols = [c for c in df.columns if c.startswith("efficiency")]

        for bus in bus_cols:
            counter = int(bus[-1])
            sdf = df.query(f"{bus} in @buses")

            if counter == 0:
                if c == "Link":
                    # add load to bus0 where energy it converted to another carrier
                    o.madd(
                        "Load",
                        sdf.index,
                        bus=sdf.bus0,
                        carrier=sdf.carrier,
                        location=sdf.location,
                        p=n.pnl(c).p0[sdf.index],
                    )
            else:
                if counter == 1:
                    p0 = n.pnl(c).p0[sdf.index]
                    p1 = n.pnl(c).p1[sdf.index]
                else:
                    efficiency = f"efficiency{counter}"
                    sdf = sdf[(sdf[efficiency] >= 0) & sdf[bus].ne("")]
                    drop = [c for c in bus_cols if not c in ["bus0", bus]] + [
                        c for c in eff_cols if not c == efficiency
                    ]
                    rename = {bus: "bus1", efficiency: "efficiency"}
                    sdf = sdf.drop(columns=drop).rename(columns=rename)
                    p0 = (
                        n.pnl(c)
                        .p0[sdf.index]
                        .rename(columns=lambda x: f"{x} {counter}")
                    )
                    p1 = n.pnl(c)[f"p{counter}"][sdf.index].rename(
                        columns=lambda x: f"{x} {counter}"
                    )
                    sdf = sdf.rename(index=lambda x: f"{x} {counter}")

                if c == "Link":
                    sdf = sdf.assign(p_nom_opt_eff=sdf.p_nom_opt * sdf.efficiency)

                if not sdf.empty:
                    o.madd(c, sdf.index, **sdf, p0=p0, p1=p1)

    comps = [
        c for c in n.one_port_components | n.branch_components if not n.df(c).empty
    ]
    carriers = list(set.union(*[set(o.df(c).carrier) for c in comps]))
    carriers = n.carriers.loc[carriers].sort_values("color")

    o.madd("Carrier", carriers.index, **carriers)
    o.kind = kind
    return o


def get_carrier_network_plotting_data(n, column):

    transmission = get_transmission_links(n)
    if column == "Optimal Capacity":
        bus_sizes = []
        bus_sizes.append(n.generators.groupby(["location", "carrier"]).p_nom_opt.sum())
        bus_sizes.append(
            n.storage_units.groupby(["location", "carrier"]).p_nom_opt.sum()
        )
        bus_sizes.append(
            n.links[~transmission].groupby(["bus1", "carrier"]).p_nom_opt_eff.sum()
        )
        bus_sizes = pd.concat(bus_sizes).fillna(0)
        bus_sizes = bus_sizes.rename(index=n.buses.location, level=0)

        line_widths = n.lines.s_nom_opt

        # only choose transmission links
        link_widths = n.links.p_nom_opt.where(transmission, 0)

        if transmission.any():
            branch_color = n.carriers.color[n.links.carrier[transmission].iloc[0]]
        else:
            branch_color = "grey"

        region_color = n.stores.groupby("location").e_nom_opt.sum()

        flow = None

    elif column == "Operation":
        bus_sizes = []

        # production
        bus_sizes.append(
            n.generators_t.p.mean()
            .groupby([n.generators.bus, n.generators.carrier])
            .sum()
        )
        bus_sizes.append(
            n.storage_units_t.p.clip(lower=0)
            .mean()
            .groupby([n.storage_units.bus, n.storage_units.carrier])
            .sum()
        )
        bus_sizes.append(
            -n.links_t.p1.mean()[~transmission]
            .groupby([n.links.bus1, n.links.carrier])
            .sum()
        )

        # consumption
        bus_sizes.append(
            -n.loads_t.p.mean().groupby([n.loads.bus, n.loads.carrier]).sum()
        )
        bus_sizes.append(
            n.storage_units_t.p.clip(upper=0)
            .mean()
            .groupby([n.storage_units.bus, n.storage_units.carrier])
            .sum()
        )
        bus_sizes.append(
            -n.links_t.p0.mean()[~transmission]
            .groupby([n.links.bus1, n.links.carrier])
            .sum()
        )

        bus_sizes = pd.concat(bus_sizes).fillna(0)
        bus_sizes = bus_sizes.rename(index=n.buses.location, level=0)

        # only choose transmission links
        flow = (
            pd.concat(
                {
                    "Line": n.lines_t.p0.mean(),
                    "Link": n.links_t.p0.mean().where(transmission, 0),
                }
            )
            * 1e2
        )
        line_widths = n.lines_t.p0.mean().abs()
        link_widths = n.links_t.p0.mean().abs().where(transmission, 0)

        if transmission.any():
            branch_color = n.carriers.color[n.links.carrier[transmission].iloc[0]]
        else:
            branch_color = "grey"

        region_color = None

    else:
        raise NotImplementedError(f"Column {column} not implemented")

    return Dict(
        bus_sizes=bus_sizes,
        line_widths=line_widths,
        link_widths=link_widths,
        region_color=region_color,
        branch_color=branch_color,
        flow=flow,
    )


def print_kind_carriers(n, carriers):
    buses = n.buses.query("carrier in @carriers").index
    for c in n.one_port_components | {"Load"}:
        df = n.df(c).query("bus in @buses")
        if df.empty:
            continue
        print(c, ":")
        for e in df.carrier.unique():
            print(" - ", e)

    for c in n.branch_components:
        print(c, ":")
        for bus in n.df(c).filter(like="bus"):
            df = n.df(c).query(f"{bus} in @buses")
            if df.empty:
                continue
            print(bus, ":")
            for e in df.carrier.unique():
                print(" - ", e)
