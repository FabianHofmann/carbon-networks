#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 13:55:01 2023

@author: fabian
"""

# Sources and sinks
# Focussing on CO2 network might be enough (turn on/off)
# How to handle import of H2?


from pathlib import Path
import warnings
import pypsa
from pypsa.descriptors import Dict
from pypsa.components import components, component_attrs
import pandas as pd
import yaml
import numpy as np

root = Path("/home/fabian/papers/co2-network")
pypsa_eur = root / "workflow/subworkflows/pypsa-eur"
pypsa_eur_sec = root / "workflow/subworkflows/pypsa-eur-sec"
override_dir = pypsa_eur_sec / "data/override_component_attrs"

SITES = [
    "residential urban decentral",
    "services urban decentral",
    "residential rural",
    "services rural",
    "services urban",
    "urban central",
]
HEATING_TECHS = ["heat", "water tanks", "solar thermal"]


CARRIER_TO_BUSES = {
    "carbon": ["co2 stored"],
    "hydrogen": ["H2"],
    "electricity": ["AC", "DC"],
    "gas": ["gas"],
}

KIND_TO_CMAP = {
    "carbon": "Reds",
    "hydrogen": "Blues",
    "electricity": "Greens",
    "gas": "Purples",
}


ONEPORT_COMPS = ["Generator", "Store", "StorageUnit"]
NODE_COMPS = ["Generator", "Link", "StorageUnit"]
REGION_COMPS = ["Store"]
BRANCH_COMPS = ["Line", "Link"]
COMPS = set(NODE_COMPS + REGION_COMPS + BRANCH_COMPS)


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

    try:
        region_cmap = KIND_TO_CMAP[n.kind]
    except AttributeError:
        region_cmap = "Blues"

    return Dict(
        bus_sizes=bus_sizes,
        line_widths=line_widths,
        link_widths=link_widths,
        region_color=region_color,
        branch_color=branch_color,
        region_cmap=region_cmap,
    )


def get_transmission_links(n):
    # only choose transmission links
    return n.links.bus0.map(n.buses.location) != n.links.bus1.map(n.buses.location)


def assign_location(n):
    for c in n.one_port_components | n.branch_components:
        df = n.df(c)

        if "location" not in df:
            df["location"] = np.nan

        ifind = pd.Series(df.index.str.find(" ", start=4), df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            if i == -1:
                continue
            names = ifind.index[ifind == i]
            df.loc[names, "location"] = names.str[:i]

        bus_col = df.columns[df.columns.str.startswith("bus")][0]
        df.location.fillna(df[bus_col], inplace=True)


def assign_interconnection(n):
    for c in n.branch_components:
        df = n.df(c)
        if df.empty:
            continue

        location0 = df.bus0.map(n.buses.location)
        location1 = df.bus1.map(n.buses.location)
        locations = pd.concat([location0, location1], axis=1)
        df["location0"] = location0
        df["location1"] = location1

        connections = locations.apply(lambda ds: " - ".join(sorted(ds)), axis=1)
        connections = connections.where(location0 != location1)
        df["interconnection"] = connections


def sanitize_locations(n):
    n.add("Bus", "EU", x=-5.5, y=46, location="EU")
    assign_location(n)
    assign_interconnection(n)
    n.buses["x"] = n.buses.location.map(n.buses.x)
    n.buses["y"] = n.buses.location.map(n.buses.y)


def fill_missing_carriers(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        new_carriers = set(c.df.carrier.unique()) - set(n.carriers.index)
        if new_carriers:
            n.madd("Carrier", list(new_carriers))


def modify_carrier_names(n):
    n.add("Carrier", "offwind", nice_name="Offshore Wind", color="#6895dd")
    n.mremove("Carrier", ["offwind-ac", "offwind-dc"])
    n.generators.loc[
        n.generators.carrier.str.startswith("offwind"), "carrier"
    ] = "offwind"
    replace = [f"{s} " for s in SITES]
    for c in n.iterate_components(
        n.one_port_components | n.branch_components | {"Bus"}
    ):
        c.df.carrier.replace(replace, "", regex=True, inplace=True)
    n.carriers.index = n.carriers.index.to_series().replace(replace, "", regex=True)
    n.carriers = n.carriers[~n.carriers.index.duplicated()]


def get_carrier_mapper(n):
    mapper = n.carriers.index.to_series()
    remove_strings = [
        " charger",
        " discharger",
        " CHP",
        " for industry",
        " transport",
        " CC",
        " to gas",
    ]
    return mapper.replace(remove_strings, "", regex=True)


def add_carrier_nice_names(n):
    nname = n.carriers.nice_name
    n.carriers.nice_name = nname.where(nname != "", nname.index.str.title())
    # replace abbreviations with capital letters
    replace = {
        "H2": "H$_2$",
        "Co2": "CO$_2$",
        "Chp": "CHP",
        "Dac": "DAC",
        "Smr": "SMR",
        " Cc": "*",
        "Ocgt": "OCGT",
        "Ac": "AC",
        "Dc": "DC",
    }
    n.carriers.nice_name.replace(replace, regex=True, inplace=True)


def add_colors(n):
    config = yaml.load(
        open(root / "config" / "config.pypsa-eur-sec.yaml"), yaml.CFullLoader
    )
    colors = config["plotting"]["tech_colors"]
    n.carriers.color = n.carriers.color.where(n.carriers.color != "")
    n.carriers.color.update(colors)
    mapper = get_carrier_mapper(n)
    n.carriers.color = n.carriers.color.fillna(mapper.map(colors))


def import_network(path):
    overrides = override_component_attrs(override_dir)
    n = pypsa.Network(path, override_component_attrs=overrides)
    sanitize_locations(n)
    fill_missing_carriers(n)
    modify_carrier_names(n)
    add_carrier_nice_names(n)
    add_colors(n)
    if n.carriers.notnull().all().all() and (n.carriers != "").all().all():
        warnings.warn("Some carriers have no color or nice_name")
    n.carriers = n.carriers.sort_values(["color"])
    return n


def override_component_attrs(directory):
    attrs = Dict({k: v.copy() for k, v in component_attrs.items()})

    for component, list_name in components.list_name.items():
        fn = directory / f"{list_name}.csv"
        if fn.exists():
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs
