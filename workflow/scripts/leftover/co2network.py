#!/usr/bin/env python
# coding: utf-8

import multiprocessing as mp
import os
import sys
from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pypsa
import seaborn as sns
import xarray as xr
import yaml
from matplotlib.legend_handler import HandlerPatch
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from pypsa.descriptors import get_switchable_as_dense as as_dense
from pypsa.plot import projected_area_factor

PATH = Path("/home/fabian/papers/co2-network/workflow/subworkflows/pypsa-eur-sec")

sys.path.append(str(PATH / "scripts"))
import holoviews as hv
from helper import override_component_attrs
from holoviews import dim, opts
from plot_network import assign_location
from plot_summary import rename_techs

hv.extension("bokeh")
hv.output(size=200)

plt.style.use(["bmh", "matplotlibrc"])
xr.set_options(display_style="html")

CLUSTERS = 181
LL = "1.0"
OPTS = "Co2L0-25H-T-H-B-I-A-solar+p3-linemaxext15-seq1000-CF+sector+co2network+false"
RUN = "20230125-carbon-small"
SCENARIO = f"elec_s_{CLUSTERS}_lv{LL}__{OPTS}_2050"
OVERRIDES = PATH / "data/override_component_attrs"

OUTPUT = "graphics-carbon"
OUTPUT_SCENARIO = f"{RUN}/{SCENARIO}/"


if not os.path.exists(OUTPUT_SCENARIO):
    os.makedirs(OUTPUT_SCENARIO)


with open(f"{RUN}/configs/config.yaml") as file:
    config = yaml.safe_load(file)


tech_colors = config["plotting"]["tech_colors"]


fn = f"resources/regions_onshore_elec_s_{CLUSTERS}.geojson"
nodes = gpd.read_file(fn).set_index("name")

fn = f"resources/regions_offshore_elec_s_{CLUSTERS}.geojson"
offnodes = gpd.read_file(fn).set_index("name")

joinednodes = pd.concat([offnodes, nodes]).dissolve(by="name")

fn = f"resources/country_shapes.geojson"
cts = gpd.read_file(fn).set_index("name")

regions = pd.concat(
    [
        gpd.read_file(f"resources/regions_onshore.geojson"),
        gpd.read_file(f"resources/regions_offshore.geojson"),
    ]
)
regions = regions.dissolve("name")

fn = f"resources/regions_onshore.geojson"
onregions = gpd.read_file(fn).set_index("name")

fn = f"resources/regions_onshore.geojson"
offregions = gpd.read_file(fn).set_index("name")

epsg = 3035
regions["Area"] = regions.to_crs(epsg=epsg).area.div(1e6)
onregions["Area"] = onregions.to_crs(epsg=epsg).area.div(1e6)
offregions["Area"] = offregions.to_crs(epsg=epsg).area.div(1e6)
nodes["Area"] = nodes.to_crs(epsg=epsg).area.div(1e6)
offnodes["Area"] = offnodes.to_crs(epsg=epsg).area.div(1e6)
joinednodes["Area"] = joinednodes.to_crs(epsg=epsg).area.div(1e6)


europe_shape = nodes.dissolve()
europe_shape.index = ["EU"]


minx, miny, maxx, maxy = europe_shape.explode(ignore_index=True).total_bounds
BOUNDARIES = [minx, maxx - 4, miny, maxy]


overrides = override_component_attrs(OVERRIDES)
fn = f"{RUN}/postnetworks/{SCENARIO}.nc"
network = pypsa.Network(fn, override_component_attrs=overrides)


unique_link_carriers = n.links.carrier.unique()
GAS_NETWORK = "gas pipeline" in unique_link_carriers
H2_NETWORK = any("H2 pipeline" in ulc for ulc in unique_link_carriers)

# ## Utilities


class HandlerCircle(HandlerPatch):
    """
    Legend Handler used to create circles for legend entries.

    This handler resizes the circles in order to match the same dimensional
    scaling as in the applied axis.
    """

    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        fig = legend.get_figure()
        ax = legend.axes

        unit = np.diff(ax.transData.transform([(0, 0), (1, 1)]), axis=0)[0][1]
        radius = orig_handle.get_radius() * unit * (72 / fig.dpi)
        center = 5 - xdescent, 3 - ydescent
        p = plt.Circle(center, radius)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


def add_legend_circles(
    ax, sizes, labels, scale=1, srid=None, patch_kw={}, legend_kw={}
):

    if srid is not None:
        area_correction = projected_area_factor(ax, n.srid) ** 2
        sizes = [s * area_correction for s in sizes]

    handles = make_legend_circles_for(sizes, scale, **patch_kw)

    legend = ax.legend(
        handles, labels, handler_map={Circle: HandlerCircle()}, **legend_kw
    )

    ax.add_artist(legend)


def add_legend_lines(ax, sizes, labels, scale=1, patch_kw={}, legend_kw={}):

    handles = [Line2D([0], [0], linewidth=s / scale, **patch_kw) for s in sizes]

    legend = ax.legend(handles, labels, **legend_kw)

    ax.add_artist(legend)


def add_legend_patch(ax, colors, labels, patch_kw={}, legend_kw={}):

    handles = [Patch(facecolor=c, **patch_kw) for c in colors]

    legend = ax.legend(handles, labels, **legend_kw)

    ax.add_artist(legend)


def make_legend_circles_for(sizes, scale=1.0, **kw):
    return [Circle((0, 0), radius=(s / scale) ** 0.5, **kw) for s in sizes]


def nodal_balance(n, carrier, time=slice(None), aggregate=None, energy=True):

    if not isinstance(carrier, list):
        carrier = [carrier]

    one_port_data = {}

    for c in n.iterate_components(n.one_port_components):

        df = c.df[c.df.bus.map(n.buses.carrier).isin(carrier)]

        if df.empty:
            continue

        s = c.pnl.p.loc[time, df.index] * df.sign

        s = s.groupby([df.bus.map(n.buses.location), df.carrier], axis=1).sum()

        one_port_data[c.list_name] = s

    branch_data = {}

    for c in n.iterate_components(n.branch_components):
        for col in c.df.columns[c.df.columns.str.startswith("bus")]:

            end = col[3:]

            df = c.df[c.df[col].map(n.buses.carrier).isin(carrier)]

            if df.empty:
                continue

            s = -c.pnl[f"p{end}"].reindex(columns=n.links.index, fill_value=0).loc[time]

            s = s.groupby([df[col].map(n.buses.location), df.carrier], axis=1).sum()

            branch_data[(c.list_name, end)] = s

    branch_balance = pd.concat(branch_data).groupby(level=[0, 2]).sum()
    one_port_balance = pd.concat(one_port_data)

    def skip_tiny(df, threshold=1e-1):
        return df.where(df.abs() > threshold)

    branch_balance = skip_tiny(branch_balance)
    one_port_balance = skip_tiny(one_port_balance)

    balance = pd.concat([one_port_balance, branch_balance]).stack(level=[0, 1])

    balance.index.set_names(["component", "bus"], level=[0, 2], inplace=True)

    if energy:
        balance = balance * n.snapshot_weightings.generators

    if aggregate is not None:
        keep_levels = balance.index.names.difference(aggregate)
        balance = balance.groupby(level=keep_levels).sum()

    return balance


# ## CO2 Network

n = network.copy()


crs = ccrs.EqualEarth()

seq = n.stores.query("carrier == 'co2 sequestered'")
joinednodes["CO2"] = seq.rename(index=seq.bus.map(n.buses.location)).e_nom_opt.div(
    1e6
)  # Mt
joinednodes["CO2"] = joinednodes["CO2"].where(joinednodes["CO2"] > 0.1)


balance = nodal_balance(n, "co2 sequestered", aggregate=["component", "snapshot"])


balance = (
    balance.where(balance > 0).drop("CO2 pipeline", level=1, errors="ignore").dropna()
)


def rename_techs_tyndp(tech):
    tech = rename_techs(tech)
    if tech == "gas":
        return "fossil gas"
    elif "oil emissions" in tech:
        return "oil emissions"
    else:
        return tech


balance = balance.unstack(0).groupby(rename_techs_tyndp).sum().T.stack().div(1e6)  # Mt


carriers = list(balance.index.levels[1])
colors = [tech_colors[c] for c in carriers]


n.madd("Carrier", carriers, color=colors)


assign_location(n)

# Drop non-electric buses so they don't clutter the plot
n.buses.drop(n.buses.index[n.buses.carrier != "AC"], inplace=True)

n.mremove("Link", n.links.query("carrier != 'CO2 pipeline'").index)

n.links.bus0 = n.links.bus0.str.replace(" co2 sequestered", "")
n.links.bus1 = n.links.bus1.str.replace(" co2 sequestered", "")


n.links["flow"] = n.snapshot_weightings.generators @ n.links_t.p0


fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": crs})

bus_size_factor = 50

n.plot(
    bus_sizes=balance / bus_size_factor,
    link_colors="k",
    branch_components=["Link"],
    ax=ax,
    geomap=True,
    link_widths=n.links.p_nom_opt.where(n.links.p_nom_opt > 1, 0) / 1000,
    flow=pd.concat({"Link": n.links.flow.where(n.links.flow.abs() > 10).div(25e4)}),
)

joinednodes = joinednodes.to_crs(crs.proj4_init)

joinednodes.plot(
    ax=ax,
    column="CO2",
    cmap="Blues",
    linewidths=0,
    legend=True,
    vmax=50,
    vmin=0,
    legend_kwds={
        "label": r"CO$_2$ sequestration [Mt/a]",
        "shrink": 0.7,
        "extend": "max",
    },
)

legend_kw = dict(
    loc="upper left",
    bbox_to_anchor=(0, 1.16),
    ncol=2,
    frameon=False,
)

add_legend_patch(ax, colors, carriers, legend_kw=legend_kw)

sizes = [10, 50]
labels = [f"{s} Mt/a" for s in sizes]
add_legend_circles(
    ax,
    sizes,
    labels,
    scale=bus_size_factor,
    srid=n.srid,
    legend_kw=dict(
        title="carbon capture", loc="upper left", bbox_to_anchor=(0.97, 1.16)
    ),
    patch_kw=dict(facecolor="lightgrey", edgecolor="k"),
)

plt.savefig(f"{OUTPUT_SCENARIO}/co2network.pdf", bbox_inches="tight")

# ## Sequestration Potential

offnodes["potential"] = (
    n.stores.e_nom_max.filter(like="co2 sequestered")
    .div(1e6)
    .rename(lambda x: x.replace(" co2 sequestered", ""))
)  # Mt


crs = ccrs.EqualEarth()

fig, ax = plt.subplots(figsize=(7, 6), subplot_kw={"projection": crs})

n.plot(
    bus_sizes=0,
    branch_components=[],
    geomap=True,
    ax=ax,
)

offnodes = offnodes.to_crs(crs.proj4_init)

offnodes.plot(
    ax=ax,
    column="potential",
    cmap="Blues",
    linewidths=0,
    legend=True,
    vmax=400,
    vmin=0,
    legend_kwds={
        "label": r"CO$_2$ sequestration potential [Mt/a]",
        "shrink": 0.7,
        "extend": "max",
    },
)

plt.savefig(f"{RUN}/sequestration_potential.pdf", bbox_inches="tight")
