import pandas as pd
import pypsa
from pypsa.descriptors import Dict
from pypsa.plot import add_legend_circles, add_legend_lines, add_legend_patches
from constants import CARRIER_TO_BUSES, ONEPORT_COMPS, BRANCH_COMPS, COMPS
from common import get_transmission_links, get_generation_carriers
from cartopy import crs as ccrs


def plot_map(
    ax,
    n,
    regions,
    column="Optimal Capacity",
    bus_scale=1.5e-5,
    branch_scale=3e-3,
    alpha=1,
    region_alpha=0.8,
    region_cmap="Blues",
):
    d = get_carrier_network_plotting_data(n, column=column)

    n.bus_scale = bus_scale
    n.branch_scale = branch_scale
    n.alpha = alpha

    # with plt.rc_context({"patch.linewidth": 0.0}):
    n.plot(
        bus_sizes=d.bus_sizes * bus_scale,
        bus_alpha=alpha,
        line_widths=d.line_widths * branch_scale,
        link_widths=d.link_widths * branch_scale,
        link_colors=d.branch_color,
        line_colors=d.branch_color,
        line_alpha=alpha,
        link_alpha=alpha,
        ax=ax,
        margin=0.2,
        color_geomap=None,
        geomap="10m",
    )
    regions = regions.assign(color=d.region_color / 1e3)
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
        linewidth=0.3,
        legend=True,
        legend_kwds={
            "label": n.kind.title() + " storage capacity [GWh]",
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
):
    legend_kwargs = {"loc": "upper left", "frameon": False}

    add_legend_circles(
        ax,
        [s * bus_scale for s in bus_sizes],
        [f"{s // 1000} GW" for s in bus_sizes],
        legend_kw={"bbox_to_anchor": (1, 1), **legend_kwargs},
    )
    add_legend_lines(
        ax,
        [s * branch_scale for s in branch_sizes],
        [f"{s // 1000} GW" for s in branch_sizes],
        legend_kw={"bbox_to_anchor": (1, 0.9), **legend_kwargs},
    )
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


def create_carrier_network(n, kind, include_eu=False):
    carriers = CARRIER_TO_BUSES[kind]
    buses = n.buses.query("carrier in @carriers").index

    o = pypsa.Network(snapshots=n.snapshots)
    o.madd("Bus", n.buses.index, **n.buses, marginal_price=n.buses_t.marginal_price)
    o.buses.loc[buses, "kind"] = True
    o.buses.kind.fillna(False, inplace=True)

    for c in ONEPORT_COMPS:
        df = n.df(c).query("bus in @buses")

        o.madd(c, df.index, **df, p=n.pnl(c).p[df.index])

    for c in BRANCH_COMPS:
        df = n.df(c)

        if not include_eu:
            df = df[df.bus0.map(n.buses.location).ne("EU")]

        bus_cols = [
            c for c in df.columns if c.startswith("bus") and not c.endswith("0")
        ]
        eff_cols = [
            c for c in df.columns if c.startswith("efficiency") and not c.endswith("0")
        ]

        for bus in bus_cols:
            counter = int(bus[-1])
            sdf = df.query(f"{bus} in @buses")

            if counter != 1:
                efficiency = f"efficiency{counter}"
                sdf = sdf[sdf[efficiency] >= 0]
                drop = [c for c in bus_cols if not c == bus] + [
                    c for c in eff_cols if not c == efficiency
                ]
                rename = {bus: "bus1", efficiency: "efficiency"}
                sdf = sdf.drop(columns=drop).rename(columns=rename)
                p0 = n.pnl(c).p0[sdf.index].rename(columns=lambda x: f"{x} {counter}")
                p1 = n.pnl(c)[f"p{counter}"][sdf.index].rename(
                    columns=lambda x: f"{x} {counter}"
                )
                sdf = sdf.rename(index=lambda x: f"{x} {counter}")
            else:
                p0 = n.pnl(c).p0[sdf.index]
                p1 = n.pnl(c).p1[sdf.index]

            if c == "Link":
                sdf = sdf.assign(p_nom_opt_eff=sdf.p_nom_opt * sdf.efficiency)
            o.madd(c, sdf.index, **sdf, p0=p0, p1=p1)

    carriers = list(set.union(*[set(o.df(c).carrier) for c in COMPS]))
    carriers = n.carriers.loc[carriers].sort_values("color")

    o.madd("Carrier", carriers.index, **carriers)
    o.kind = kind
    return o


def get_carrier_network_plotting_data(n, column):

    transmission = get_transmission_links(n)
    if column == "Optimal Capacity":
        bus_sizes = []
        bus_sizes.append(n.generators.groupby(["bus", "carrier"]).p_nom_opt.sum())
        bus_sizes.append(n.storage_units.groupby(["bus", "carrier"]).p_nom_opt.sum())
        bus_sizes.append(
            n.links[~transmission].groupby(["bus1", "carrier"]).p_nom_opt_eff.sum()
        )
        bus_sizes = pd.concat(bus_sizes).fillna(0)

        line_widths = n.lines.s_nom_opt

        # only choose transmission links
        link_widths = n.links.p_nom_opt.where(transmission, 0)

        branch_color = n.carriers.color[
            n.links.carrier[get_transmission_links(n)].iloc[0]
        ]

        region_color = n.stores.groupby("location").e_nom_opt.sum()
    else:
        raise NotImplementedError(f"Column {column} not implemented")

    return Dict(
        bus_sizes=bus_sizes,
        line_widths=line_widths,
        link_widths=link_widths,
        region_color=region_color,
        branch_color=branch_color,
    )
