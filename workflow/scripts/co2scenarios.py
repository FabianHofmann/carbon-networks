#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys
from itertools import product
from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import IndexSlice as idx
import xarray as xr
import yaml
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from matplotlib.transforms import Bbox
from vresutils.costdata import annuity


PATH = Path("/home/fabian/papers/co2-network/workflow/subworkflows/pypsa-eur-sec")

sys.path.append(
    "/home/fabian/papers/co2-network/workflow/subworkflows/pypsa-eur-sec/scripts"
)
from plot_summary import rename_techs

CLUSTERS = 181
OUTPUT = "20230125-carbon-small/"

MAIN_SCENARIOS = Path("20230125-carbon-small")

plt.style.use(["bmh", "matplotlibrc"])

with open(MAIN_SCENARIOS / "configs/config.yaml") as file:
    config = yaml.safe_load(file)

if not os.path.exists(OUTPUT):
    os.makedirs(OUTPUT)


def rename_techs_tyndp(tech):
    tech = rename_techs(tech)
    if "heat pump" in tech or "resistive heater" in tech:
        return "power-to-heat"
    elif tech in ["H2 Electrolysis"]:  # , "H2 liquefaction"]:
        return "power-to-hydrogen"
    elif "H2 pipeline" in tech:
        return "H2 pipeline"
    elif tech == "H2":
        return "H2 storage"
    elif tech in ["OCGT", "CHP", "gas boiler", "H2 Fuel Cell"]:
        return "gas-to-power/heat"
    # elif "solar" in tech:
    #    return "solar"
    elif tech in ["Fischer-Tropsch", "methanolisation"]:
        return "power-to-liquid"
    elif "offshore wind" in tech:
        return "offshore wind"
    elif "SMR" in tech:
        return tech.replace("SMR", "steam methane reforming")
    elif "DAC" in tech:
        return "direct air capture"
    elif "CC" in tech or "sequestration" in tech:
        return "carbon capture"
    else:
        return tech


preferred_order = pd.Index(
    [
        "transmission lines",
        "electricity distribution grid",
        "fossil oil and gas",
        "hydroelectricity",
        "hydro reservoir",
        "run of river",
        "pumped hydro storage",
        "solid biomass",
        "biogas",
        "onshore wind",
        "offshore wind",
        "offshore wind (AC)",
        "offshore wind (DC)",
        "solar PV",
        "solar thermal",
        "solar rooftop",
        "solar",
        "building retrofitting",
        "ground heat pump",
        "air heat pump",
        "heat pump",
        "resistive heater",
        "power-to-heat",
        "gas-to-power/heat",
        "CHP",
        "OCGT",
        "gas boiler",
        "gas",
        "natural gas",
        "helmeth",
        "methanation",
        "power-to-gas",
        "power-to-hydrogen",
        "H2 pipeline",
        "H2 liquefaction",
        "H2 storage",
        "hydrogen storage",
        "power-to-liquid",
        "battery storage",
        "hot water storage",
        "CO2 sequestration",
        "CCS",
        "carbon capture and sequestration",
        "DAC",
        "direct air capture",
    ]
)


def rename_techs_balances(tech):
    tech = rename_techs(tech)
    if "heat pump" in tech:
        return "heat pump"
    elif tech in ["H2 Electrolysis"]:  # , "H2 liquefaction"]:
        return "power-to-hydrogen"
    elif "solar" in tech:
        return "solar"
    elif tech in ["Fischer-Tropsch", "methanolisation"]:
        return "power-to-liquid"
    elif tech == "DAC":
        return "direct air capture"
    elif "offshore wind" in tech:
        return "offshore wind"
    elif tech == "oil" or tech == "gas":
        return "fossil oil and gas"
    elif tech in ["BEV charger", "V2G", "Li ion", "land transport EV"]:
        return "battery electric vehicles"
    elif tech in ["biogas", "solid biomass"]:
        return "biomass"
    elif tech in ["industry electricity", "electricity", "agriculture electricity"]:
        return "electricity demand"
    elif tech in ["agriculture heat", "heat", "low-temperature heat for industry"]:
        return "heat demand"
    elif "solid biomass for industry" in tech:
        return "biomass demand"
    elif "gas for industry" in tech:
        return "methane demand"
    elif tech in ["H2 for industry", "land transport fuel cell"]:
        return "hydrogen demand"
    elif tech in [
        "kerosene for aviation",
        "naphtha for industry",
        "shipping methanol",
        "agriculture machinery oil",
    ]:
        return "liquid hydrocarbon demand"
    elif tech in [
        "transmission lines",
        "H2 pipeline",
        "H2 pipeline retrofitted",
        "H2",
        "electricity distribution grid",
        "hot water storage",
        "SMR",
        "SMR CC",
        "OCGT",
        "CHP",
        "gas boiler",
        "H2 Fuel Cell",
        "resistive heater",
        "battery storage",
        "methanation",
    ]:
        return "other"
    else:
        return tech


def parse_index(c, with_resolution=False):

    match = re.search(r"seq([0-9.]*)", c[2])
    seq = 200 if match is None else float(match.groups()[0])

    h2 = "no H2 grid" if "noH2network" in c[2] else "H2 grid"

    co2 = "no CO2 grid" if "co2network+false" in c[2] else "CO2 grid"

    to_return = (seq, h2, co2)

    return to_return


def load_main(scenarios=None, rename=True):

    if scenarios is None:
        scenarios = MAIN_SCENARIOS

    costs = pd.read_csv(
        scenarios / "csvs/costs.csv", header=[0, 1, 2, 3], index_col=[0, 1, 2]
    )

    costs = costs.xs(str(2050), level="planning_horizon", axis=1)

    names = ["seq", "h2", "co2"]
    costs.columns = pd.MultiIndex.from_tuples(
        [parse_index(c) for c in costs.columns], names=names
    )

    df = costs.groupby(level=2).sum().div(1e9)

    if rename:
        df = df.groupby(df.index.map(rename_techs_tyndp)).sum()

    to_drop = df.index[df.max(axis=1).fillna(0.0) < 1.2]
    print(to_drop)
    df.drop(to_drop, inplace=True)

    order = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )
    df = df.loc[order]

    return df


def load_main_supply_energy(
    scenarios=None,
    rename=True,
    carrier="energy",
):

    if scenarios is None:
        scenarios = MAIN_SCENARIOS

    df = pd.read_csv(
        scenarios / "csvs/supply_energy.csv", index_col=[0, 1, 2], header=[0, 1, 2, 3]
    )

    co2_carriers = ["co2", "co2 stored", "process emissions"]
    if carrier == "energy":
        carrier = [i for i in df.index.levels[0] if i not in co2_carriers]

    df = df.loc[carrier].groupby(level=2).sum().div(1e6)  # TWh / MtCO2
    df.index = [
        i[:-1]
        if ((i not in ["co2", "NH3", "H2"]) and (i[-1:] in ["0", "1", "2", "3"]))
        else i
        for i in df.index
    ]

    df = df.xs("2050", level="planning_horizon", axis=1)

    names = ["seq", "h2", "co2"]

    df.columns = pd.MultiIndex.from_tuples(
        [parse_index(c) for c in df.columns], names=names
    )

    if rename:
        df = df.groupby(df.index.map(rename_techs_tyndp)).sum()

    to_drop = df.index[df.abs().max(axis=1).fillna(0.0) < 10]
    df.drop(to_drop, inplace=True)

    order = preferred_order.intersection(df.index).append(
        df.index.difference(preferred_order)
    )
    df = df.loc[order]

    return df


def add_arrow(ax0, ax1, pos0, pos1, **arrow_kwargs):
    ax0tr = ax0.transData  # Axis 0 -> Display
    ax1tr = ax1.transData  # Axis 1 -> Display
    figtr = fig.transFigure.inverted()  # Display -> Figure
    ptB = figtr.transform(ax0tr.transform(pos0))
    ptE = figtr.transform(ax1tr.transform(pos1))
    arrow = FancyArrowPatch(
        ptB,
        ptE,
        transform=fig.transFigure,  # Place arrow in figure coord system
        **arrow_kwargs,
    )
    fig.patches.append(arrow)


def rename_techs_carbon_balances(tech):
    prefix_to_remove = [
        "residential ",
        "services ",
        "urban ",
        "rural ",
        "central ",
        "decentral ",
    ]
    for ptr in prefix_to_remove:
        if tech[: len(ptr)] == ptr:
            tech = tech[len(ptr) :]
    if tech in [
        "oil emissions",
        "agriculture machinery oil emissions",
        "shipping methanol emissions",
    ]:
        return "liquid hydrocarbons emissions"
    elif tech == "biogas to gas":
        return "biogas upgrading"
    # elif tech == "DAC":
    #    return "direct air capture"
    elif "SMR" in tech:
        return tech.replace("SMR", "steam methane reforming")
    else:
        return tech


tech_colors = config["plotting"]["tech_colors"]
tech_colors = pd.Series(tech_colors).rename(rename_techs_tyndp).to_dict()
tech_colors["battery electric vehicles"] = tech_colors["BEV charger"]
tech_colors["other"] = "#454545"
tech_colors["heat demand"] = "indianred"
tech_colors["electricity demand"] = tech_colors["electricity"]
tech_colors["hydrogen demand"] = tech_colors["land transport fuel cell"]
tech_colors["methane demand"] = tech_colors["helmeth"]
tech_colors["liquid hydrocarbon demand"] = tech_colors["kerosene for aviation"]
tech_colors["biomass demand"] = tech_colors["biogas"]
tech_colors["fossil oil and gas"] = tech_colors["fossil gas"]
tech_colors["biogas upgrading"] = tech_colors["biogas"]
tech_colors["liquid hydrocarbons emissions"] = tech_colors["kerosene for aviation"]
tech_colors["solar"] = tech_colors["solar PV"]
tech_colors["heat pump"] = "#2fb537"


# ## H2 Network Scenarios

SCENARIO = MAIN_SCENARIOS

df = load_main(SCENARIO)

df.sum().sort_values()

df.sum().unstack("co2")

seq = 200
df_new = df.xs(seq, level="seq", axis=1)

dff = df_new.sum().unstack("h2")

dff / dff.min(axis=0)

dff.sort_index(ascending=False, inplace=True)

h2_rel_benefit = (dff.T / dff.min(axis=1) * 100 - 100).iloc[1].reset_index(drop=True)
h2_abs_benefit = (dff.T - dff.min(axis=1)).iloc[1].reset_index(drop=True)

ac_rel_benefit = (dff / dff.min(axis=0) * 100 - 100).iloc[0].reset_index(drop=True)
ac_abs_benefit = (dff - dff.min(axis=0)).iloc[0].reset_index(drop=True)

max_rel_benefit = dff.max().max() / dff.min().min() * 100 - 100

max_abs_benefit = int(dff.max().max() - dff.min().min())


xx = enumerate(df_new.columns.get_level_values("co2").unique())  # [::-1])
yy = enumerate(df_new.columns.get_level_values("h2").unique())

fig, axs = plt.subplots(2, 2, figsize=(4, 6), sharey=True)

plt.subplots_adjust(hspace=0.5, wspace=1)

kwargs = dict(stacked=True, color=tech_colors, ylim=(0, 900), legend=False)

for x, y in product(xx, yy):

    ax = axs[x[0], y[0]]
    toprow_kwargs = (
        dict(
            title="with\nhydrogen grid\n"
            if y[1] == "H2 grid"
            else "without\nhydrogen grid\n"
        )
        if x[0] == 0
        else {}
    )
    ylabel = (
        "with\ncarbon-dioxide grid\n\nbn€/a"
        if x[1] == "CO2 grid"
        else "without\ncarbon-dioxide grid\n\nbn€/a"
    )

    to_plot = df_new.xs((x[1], y[1]), axis=1, level=["co2", "h2"]).T
    to_plot.plot.bar(
        ax=ax,
        ylabel=ylabel,
        **kwargs,
        **toprow_kwargs,
    )
    ax.set_xlabel("", rotation=0)
    ax.set_xticks([], [])
    ax.tick_params(labelrotation=0)
    ax.grid(axis="y")
    ax.title.set_size(11)
    ax.set_yticks(np.arange(0, 901, 100))
    ax.text(-0.3, 825, f"{to_plot.sum().sum():.0f} bn€/a", color="grey", fontsize=9.5)

    for i in ["top", "right", "left", "bottom"]:
        ax.spines[i].set_visible(False)

handles, labels = ax.get_legend_handles_labels()
handles.reverse()
labels.reverse()
fig.legend(handles, labels, bbox_to_anchor=(1.55, 1))

fig.text(
    0.41,
    0.3,
    f"+ {h2_rel_benefit[0]:.1f}%\n+ {h2_abs_benefit[0]:.0f} bn€/a",
    fontsize=11,
)
fig.text(
    0.41,
    0.76,
    f"+ {h2_rel_benefit[1]:.1f}%\n+ {h2_abs_benefit[1]:.0f} bn€/a",
    fontsize=11,
)
fig.text(
    0.11,
    0.47,
    f"+ {ac_rel_benefit[0]:.1f}%\n+ {ac_abs_benefit[0]:.0f}\nbn€/a",
    fontsize=11,
)
fig.text(
    0.8,
    0.47,
    f"+ {ac_rel_benefit[1]:.1f}%\n+ {ac_abs_benefit[1]:.0f}\nbn€/a",
    fontsize=11,
)
fig.text(
    0.33,
    0.47,
    f"+ {max_rel_benefit:.1f}%\n+ {max_abs_benefit} bn€/a",
    fontsize=11,
    color="grey",
)


norm = mpl.colors.Normalize(vmin=0, vmax=10)
m = cm.ScalarMappable(norm=norm, cmap=cm.cividis)

arrow_style = dict(arrowstyle="simple", mutation_scale=22, ec=None)

add_arrow(
    axs[0, 0],
    axs[0, 1],
    (0.5, 500),
    (-0.5, 500),
    fc=m.to_rgba(h2_rel_benefit[1]),
    **arrow_style,
)

add_arrow(
    axs[1, 0],
    axs[1, 1],
    (0.5, 500),
    (-0.5, 500),
    fc=m.to_rgba(h2_rel_benefit[0]),
    **arrow_style,
)

add_arrow(
    axs[0, 0],
    axs[1, 0],
    (0, 0),
    (0, 900),
    fc=m.to_rgba(ac_rel_benefit[0]),
    **arrow_style,
)

add_arrow(
    axs[0, 1],
    axs[1, 1],
    (0, 0),
    (0, 900),
    fc=m.to_rgba(ac_rel_benefit[1]),
    **arrow_style,
)

add_arrow(
    axs[0, 0],
    axs[1, 1],
    (0.5, 0),
    (-0.5, 900),
    fc=m.to_rgba(max_rel_benefit),
    **arrow_style,
)

plt.savefig(f"{OUTPUT}sensitivity-h2-co2-seq{seq}.pdf", bbox_inches="tight")


# %% cost plot

df = load_main()
df = df.xs("H2 grid", level="h2", axis=1)[[200, 1000]].T
df.index = ["{}Mt seq. \n {}".format(int(i), j) for i, j in df.index]

fig, ax = plt.subplots(figsize=(6, 4))
colors = [tech_colors[i] for i in df.columns]
df.plot.bar(ax=ax, stacked=True, linewidth=0, color=colors)

handles, labels = ax.get_legend_handles_labels()
handles.reverse()
labels.reverse()

# ax.set_xlim(200, 1000)
# ax.set_xlabel(r"Sequestration Potential [Mt$_{CO_2}$/a]")
ax.set_ylabel("System Cost [bn€/a]")
ax.grid(axis="y")
# legend on side
ax.legend(
    handles, labels, ncol=2, frameon=False, bbox_to_anchor=(0, 1.05), loc="lower left"
)

for i in ["top", "right", "left", "bottom"]:
    ax.spines[i].set_visible(False)

# ax.set_axisbelow(False)
# ax.tick_params(rotation=0, labelsize=10)
# ax.set_xticks([200, 400, 600, 800, 1000])
# fig.tight_layout()
fig.savefig(OUTPUT + "seq-sensitivity-co2-h2.pdf", bbox_inches="tight")

# %%

cols = ["no H2 grid", "no CO2 grid"]

scenarios = MAIN_SCENARIOS
co2_carriers = ["co2", "co2 stored", "process emissions"]
balances_df = pd.read_csv(
    scenarios / "csvs/supply_energy.csv", index_col=[0, 1, 2], header=[0, 1, 2, 3]
)

balances = {i.replace(" ", "_"): [i] for i in balances_df.index.levels[0]}
balances["energy"] = [i for i in balances_df.index.levels[0] if i not in co2_carriers]
balances["carbon"] = [i for i in balances_df.index.levels[0] if i in co2_carriers]

key = "energy"

df = balances_df.loc[balances[key]]
df = df.groupby(level=2).sum().div(1e6)
df.index = [
    i[:-1]
    if ((i not in ["co2", "NH3", "H2"]) and (i[-1:] in ["0", "1", "2", "3"]))
    else i
    for i in df.index
]


df = df.groupby(rename_techs_balances).sum()
df.columns = pd.MultiIndex.from_tuples(
    [parse_index(c) for c in df.columns], names=["seq", "h2", "co2"]
)

df = df.loc[:, idx[[200, 1000], "H2 grid"]].droplevel(1, axis=1)

order = pd.Index(
    [
        "electricity demand",
        "battery electric vehicles",
        "heat demand",
        "hydrogen demand",
        "biomass demand",
        "methane demand",
        "liquid hydrocarbon demand",
        "power-to-liquid",
        "methanation",
        "power-to-hydrogen",
        "other",
        "direct air capture",
        "fossil oil and gas",
        "hydroelectricity",
        "biomass",
        "offshore wind",
        "onshore wind",
        "solar",
        "heat pump",
    ]
)

order = order.intersection(df.index).append(df.index.difference(order))
df = df.loc[order]
e = df.copy()


# ## Carbon Balance

key = "co2"
df = balances_df.loc[balances[key]]
df = df.groupby(level=2).sum().div(1e6)
df.index = [
    i[:-1]
    if ((i not in ["co2", "NH3", "H2"]) and (i[-1:] in ["0", "1", "2", "3"]))
    else i
    for i in df.index
]
df = df.groupby(rename_techs_carbon_balances).sum()
df.columns = pd.MultiIndex.from_tuples(
    [parse_index(c) for c in df.columns], names=["seq", "h2", "co2"]
)
df = df.loc[:, idx[[200, 1000], "H2 grid"]].droplevel(1, axis=1)
df.drop("co2", inplace=True)

order = pd.Index(
    [
        "liquid hydrocarbons emissions",
        "methanol emissions",
        "process emissions CC",
        "gas for industry CC",
        "gas CHP CC",
        "gas CHP",
        "OCGT",
        "gas boiler",
        "steam methane reforming",
        "steam methane reforming CC",
        "biogas upgrading",
        "solid biomass CHP CC",
        "solid biomass for industry CC",
        "DAC",
    ]
)

order = order.intersection(df.index).append(df.index.difference(order))
df = df.loc[order]

df = df.loc[df.abs().max(axis=1) > 0.01]

tech_colors["biogas upgrading"] = tech_colors["biogas"]
tech_colors["liquid hydrocarbons emissions"] = tech_colors["kerosene for aviation"]

df = df.sort_index(axis=1, level=1)  # .droplevel(1, axis=1)
e = e.sort_index(axis=1, level=1)  # .droplevel(1, axis=1)


e.columns = ["{}Mt seq. \n {}".format(int(i), j) for i, j in e.columns]
df.columns = ["{}Mt seq. \n {}".format(int(i), j) for i, j in df.columns]


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 4.75))

lim = 15000
e.T.plot.bar(ax=ax1, stacked=True, color=tech_colors, ylim=(-lim, lim))

lim = 800
df.T.plot.bar(ax=ax2, stacked=True, cmap="tab20", ylim=(-lim, lim))

ax1.legend(loc=(1.03, -0.1), ncol=1)
ax2.legend(loc=(1.03, -0.0))

ax1.grid(axis="y")
ax2.grid(axis="y")

ax1.set_xlabel(None)
ax2.set_xlabel(None)

ax1.set_ylabel(r"consumption $\leftarrow$ TWh $\rightarrow$ supply             ")
ax2.set_ylabel(r"withdrawal $\leftarrow$ Mt$_{CO_2}$ $\rightarrow$ emission      ")

# ax1.set_xlabel("sequestration potential [Mt$_{CO_2}$/a]")
# ax2.set_xlabel("sequestration potential [Mt$_{CO_2}$/a]")

ax1.set_title("Energy Balance", fontsize=11)
ax2.set_title(r"CO$_2$ Balance (atmosphere)", fontsize=11)

plt.tight_layout()

plt.savefig(OUTPUT + "balance.pdf", bbox_inches="tight")
