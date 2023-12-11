#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 13:55:01 2023

@author: fabian
"""
import os
from pathlib import Path
import textwrap
import warnings
import pypsa
from pypsa.descriptors import Dict
from pypsa.components import components, component_attrs
import pandas as pd
from snakemake.utils import update_config
import yaml
import numpy as np
import logging

logger = logging.getLogger(__name__)

root = Path(__file__).parent.parent.parent.resolve()
pypsa_eur = root / "workflow/subworkflows/pypsa-eur"


SITES = [
    "residential urban decentral",
    "services urban decentral",
    "residential rural",
    "services rural",
    "services urban",
    "urban central",
]


def set_scenario_config(snakemake):
    if snakemake.config["run"]["scenarios"] and "run" in snakemake.wildcards:
        try:
            with open(snakemake.config["scenariofile"], "r") as f:
                scenario_config = yaml.safe_load(f)
        except FileNotFoundError:
            # fallback for mock_snakemake
            script_dir = Path(__file__).parent.resolve()
            root_dir = script_dir.parent
            with open(root_dir / snakemake.config["scenariofile"], "r") as f:
                scenario_config = yaml.safe_load(f)
        update_config(snakemake.config, scenario_config[snakemake.wildcards.run])


def import_network(
    path,
    revert_dac=False,
    offshore_sequestration=False,
    offshore_regions=None,
    remove_gas_store_capital_cost=False,
):
    n = pypsa.Network(path)
    if offshore_sequestration:
        move_sequestration_offshore(n, offshore_regions)
    if revert_dac:
        revert_dac_ports(n)
    if remove_gas_store_capital_cost:
        remove_gas_store_capex(n)
    sanitize_locations(n)
    fill_missing_carriers(n)
    modify_carrier_names(n)
    add_carrier_nice_names(n)
    add_carrier_groups(n)
    update_colors(n)
    n.carriers.drop("", inplace=True)
    if (no_names := n.carriers.query("nice_name == ''").index).any():
        warnings.warn(f"Carriers {no_names} have no nice_name")
    if (no_colors := n.carriers.query("color == ''").index).any():
        warnings.warn(f"Carriers {no_colors} have no color")
    if (no_groups := n.carriers.query("group == ''").index).any():
        warnings.warn(f"Carriers {no_groups} have no technology group")
    if (no_group_colors := n.carriers.query("group_color == ''").index).any():
        warnings.warn(f"Carriers {no_group_colors} have no technology group color")
    n.carriers = n.carriers.sort_values(["color"])
    return n


def revert_dac_ports(n):
    dac = n.links.index[n.links.carrier == "DAC"]
    n.links.loc[dac, ["bus0", "bus2"]] = n.links.loc[dac, ["bus2", "bus0"]].values
    n.links_t.p0[dac], n.links_t.p2[dac] = (
        n.links_t.p2[dac],
        n.links_t.p0[dac].values,
    )


def move_sequestration_offshore(n, offshore_regions):
    if offshore_regions is None:
        raise ValueError(
            "offshore regions must be given for reallocating sequestration"
        )

    offcoords = pd.concat(
        [offshore_regions.centroid.x, offshore_regions.centroid.y],
        axis=1,
        keys=["x", "y"],
    )

    hubs = offcoords.index
    hubs_co2 = hubs + " co2 stored"

    offbuses = n.madd(
        "Bus",
        offcoords.index,
        suffix=" offshore",
        location=offcoords.index + " offshore",
        carrier="AC",
        **offcoords,
    )
    offbuses_co2 = n.madd(
        "Bus",
        offcoords.index,
        suffix=" co2 stored offshore",
        carrier="co2 stored",
        location=offcoords.index + " co2 stored offshore",
        **offcoords,
    )

    sequestration_buses = n.buses.query("carrier == 'co2 sequestered'").index
    n.buses.loc[sequestration_buses, "location"] += " offshore"
    sequestration_links = n.links.query("carrier == 'co2 sequestered'").index
    n.links.loc[sequestration_links, "bus0"] += " offshore"

    rename = lambda col: col.replace(" co2 sequestered", "")
    p = n.links_t.p0[sequestration_links].rename(columns=rename)[hubs].values
    names = "CO2 pipeline " + hubs + " -> " + offbuses
    n.madd(
        "Link",
        names,
        bus0=hubs_co2,
        bus1=offbuses_co2,
        carrier="CO2 pipeline",
        p0=p,
        p1=-p,
    )


def remove_gas_store_capex(n):
    """Remove capex contributions from gas storages"""
    gas_stores = n.stores.filter(like="gas Store", axis=0).index
    n.stores.loc[gas_stores, "capital_cost"] = 0


def group_small_contributions(df, threshold=0.01):
    if isinstance(df, pd.DataFrame):
        is_small = df.sum() < threshold
        if is_small.sum() <= 1:
            return df
        small = df.loc[:, is_small]
        df = df.loc[:, ~is_small]
        df = df.assign(Other=small.sum(1))
        return df
    elif isinstance(df, pd.Series):
        is_small = df < threshold
        if is_small.sum() <= 1:
            return df
        small = df[is_small]
        df = df[~is_small]
        df.assign(Other=small.sum())
        return df
    else:
        raise TypeError("df must be a DataFrame or Series")


def sort_rows_by_diff(df: pd.DataFrame):
    means = df.mean(axis=1)
    df_pos = df[means > 0]
    df_neg = df[means <= 0]

    col_pos = df_pos.max(axis=1) - df_pos.min(axis=1)
    col_neg = df_neg.max(axis=1) - df_neg.min(axis=1)

    df_pos = (
        df_pos.assign(diff=col_pos)
        .sort_values(by="diff", ascending=True)
        .drop("diff", axis=1)
    )
    df_neg = (
        df_neg.assign(diff=col_neg)
        .sort_values(by="diff", ascending=True)
        .drop("diff", axis=1)
    )

    return pd.concat([df_neg, df_pos])


def sort_rows_by_relative_diff(df: pd.DataFrame):
    means = df.mean(axis=1)
    df_pos = df[means > 0]
    df_neg = df[means <= 0]

    col_pos = (df_pos.max(axis=1) - df_pos.min(axis=1)) / df_pos.mean(axis=1)
    col_neg = (df_neg.max(axis=1) - df_neg.min(axis=1)) / df_neg.mean(axis=1).abs()

    df_pos = (
        df_pos.assign(diff=col_pos)
        .sort_values(by="diff", ascending=True)
        .drop("diff", axis=1)
    )
    df_neg = (
        df_neg.assign(diff=col_neg)
        .sort_values(by="diff", ascending=True)
        .drop("diff", axis=1)
    )

    return pd.concat([df_neg, df_pos])


def groupby_carrier_across_cc(n, c, nice_names) -> pd.Series:
    return n.df(c).carrier.replace(" CC", "", regex=True).replace(n.carriers.nice_name)


def get_ordered_handles_labels(ax, data, wrap=20):
    """
    Use this for stacked bar plots with negative values.
    """
    handles, labels = ax.get_legend_handles_labels()
    legend_data = pd.Series(handles, index=labels)
    is_negative = data.mean(1).lt(0)
    res = pd.concat([legend_data[~is_negative][::-1], legend_data[is_negative]])
    res = res[~res.index.duplicated()]
    if wrap:
        wrapper = textwrap.TextWrapper(width=wrap)
        res.rename(wrapper.fill, inplace=True)

    return res.values, res.index


def assert_carriers_existent(n, carriers, c):
    if c == "Link":
        if not n.meta["sector"]["co2network"]:
            carriers = set(carriers) - {"CO2 pipeline"}
        if not n.meta["sector"]["H2_network"]:
            carriers = set(carriers) - {"H2 pipeline"}
        if not n.meta["sector"]["gas_network"]:
            carriers = set(carriers) - {"gas pipeline", "gas pipeline new"}
    if not set(carriers).issubset(n.df(c).carrier.unique()):
        logger.warning(
            f"Carriers {set(carriers) - set(n.df(c).carrier)} are not in the component {c}."
        )


def get_transmission_links(n, with_eu=False):
    # only choose transmission links
    if with_eu:
        return n.links.bus0.map(n.buses.location) != n.links.bus1.map(n.buses.location)
    return (
        (n.links.bus0.map(n.buses.location) != n.links.bus1.map(n.buses.location))
        & ~n.links.bus0.map(n.buses.location).str.contains("EU")
        & ~n.links.bus1.map(n.buses.location).str.contains("EU")
    )


def get_p_nom_opt_eff(links, bus):
    port = int(bus[-1])
    p_nom_opt = links.p_nom_opt
    if port == 0:
        efficiency = 1
    elif port == 1:
        efficiency = links.efficiency.abs()
    else:
        efficiency = links["efficiency" + str(port)].abs()
    return p_nom_opt * efficiency


def get_carrier_production(n, kind, config, which="capacity"):
    """
    Get the production of a certain kind of carrier.

    The corresponding values are calculated on the basis of the config file entry
    'kind_to_carrier'.

    Parameters
    ----------
    n : pypsa.Network
    kind : str
    config : dict
    which: str, optional
        Whether to return the capacity or the operation of the carriers, by default 'capacity'
    """
    weights = n.snapshot_weightings.generators
    location = n.buses.location
    specs = config["constants"]["kind_to_carrier"][kind].get("production", {})
    res = []

    carriers = specs.get("Generator", [])
    assert_carriers_existent(n, carriers, "Generator")
    gens = n.generators.query("carrier in @carriers")
    groups = [gens.bus.map(n.buses.location), gens.carrier]
    if which == "capacity":
        df = gens.groupby(groups).p_nom_opt.sum()
    elif which == "operation":
        df = (weights @ n.generators_t.p)[gens.index].groupby(groups).sum()
    res.append(df)

    buses = specs.get("Link", [])
    for bus, carriers in buses.items():
        port = int(bus[-1])
        links = n.links.query("carrier in @carriers")
        assert_carriers_existent(n, carriers, "Link")
        groups = [links[bus].map(location), links.carrier]
        if which == "capacity":
            df = get_p_nom_opt_eff(links, bus).groupby(groups).sum()
        elif which == "operation":
            p = f"p{port}"
            # producing operation at a port is always negative
            df = (
                -(weights @ n.links_t[p][links.index].clip(upper=0))
                .groupby(groups)
                .sum()
            )
        res.append(df)

    carriers = specs.get("StorageUnit", [])
    assert_carriers_existent(n, carriers, "StorageUnit")
    stos = n.storage_units.query("carrier in @carriers")
    groups = [stos.bus.map(n.buses.location), stos.carrier]
    if which == "capacity":
        df = stos.groupby(groups).p_nom_opt.sum()
    elif which == "operation":
        df = (
            (weights @ n.storage_units_t.p.clip(lower=0))[stos.index]
            .groupby(groups)
            .sum()
        )
    res.append(df)

    carriers = specs.get("Store", [])
    assert_carriers_existent(n, carriers, "Store")
    sus = n.stores.query("carrier in @carriers")
    groups = [sus.bus.map(location), sus.carrier]
    if which == "operation":
        df = (weights @ n.stores_t.p[sus.index].clip(lower=0)).groupby(groups).sum()
        res.append(df)

    carriers = specs.get("Load", [])
    assert_carriers_existent(n, carriers, "Load")
    loads = n.loads.query("carrier in @carriers")
    groups = [loads.bus.map(n.buses.location), loads.carrier]
    if which == "operation":
        p_set = n.get_switchable_as_dense("Load", "p_set")
        df = -(weights @ p_set[loads.index].clip(upper=0)).groupby(groups).sum()
    res.append(df)

    return pd.concat(res)


def get_carrier_consumption(n, kind, config, which="capacity"):
    """
    Get the consumption of a certain kind of carrier.

    The corresponding values are calculated on the basis of the config file entry
    'kind_to_carrier'.

    Parameters
    ----------
    n : pypsa.Network
    kind : str
    config : dict
    which: str, optional
        Whether to return the capacity or the operation of the carriers, by default 'capacity'
    """
    weights = n.snapshot_weightings.generators
    location = n.buses.location
    specs = config["constants"]["kind_to_carrier"][kind].get("consumption", {})
    res = []

    carriers = specs.get("Load", [])
    assert_carriers_existent(n, carriers, "Load")
    loads = n.loads.query("carrier in @carriers")
    groups = [loads.bus.map(location), loads.carrier]
    if which == "operation":
        p_set = n.get_switchable_as_dense("Load", "p_set")
        df = (weights @ p_set[loads.index].clip(lower=0)).groupby(groups).sum()
        res.append(df)

    buses = specs.get("Link", [])
    for bus, carriers in buses.items():
        port = int(bus[-1])
        assert_carriers_existent(n, carriers, "Link")
        links = n.links.query("carrier in @carriers")
        groups = [links[bus].map(location), links.carrier]
        if which == "capacity":
            df = get_p_nom_opt_eff(links, bus).groupby(groups).sum()
        elif which == "operation":
            p = f"p{port}"
            df = (
                (weights @ n.links_t[p][links.index].clip(lower=0))
                .groupby(groups)
                .sum()
            )
        res.append(df)

    carriers = specs.get("StorageUnit", [])
    assert_carriers_existent(n, carriers, "StorageUnit")
    stos = n.storage_units.query("carrier in @carriers")
    groups = [stos.bus.map(location), stos.carrier]
    if which == "capacity":
        df = stos.groupby(groups).p_nom_opt.sum()
    elif which == "operation":
        df = (
            -(weights @ n.storage_units_t.p[stos.index].clip(upper=0))
            .groupby(groups)
            .sum()
        )
    res.append(df)

    carriers = specs.get("Store", [])
    assert_carriers_existent(n, carriers, "Store")
    sus = n.stores.query("carrier in @carriers")
    groups = [sus.bus.map(location), sus.carrier]
    if which == "capacity":
        df = sus.groupby(groups).e_nom_opt.sum()
    elif which == "operation":
        df = -(weights @ n.stores_t.p[sus.index].clip(upper=0)).groupby(groups).sum()
    res.append(df)

    return pd.concat(res)


def get_carrier_storage(n, kind, config, which="capacity"):
    """
    Get the storage of a certain kind of carrier.

    The corresponding values are calculated on the basis of the config file entry
    'kind_to_carrier'.

    Parameters
    ----------
    n : pypsa.Network
    kind : str
    config : dict
    which: str, optional
        Whether to return the capacity or the operation of the carriers, by default 'capacity'
    """
    weights = n.snapshot_weightings.generators
    specs = config["constants"]["kind_to_carrier"][kind].get("storage", {})
    res = []

    location = n.buses.location
    carriers = specs.get("StorageUnit", [])
    assert_carriers_existent(n, carriers, "StorageUnit")
    stos = n.storage_units.query("carrier in @carriers")
    groups = [stos.bus.map(location), stos.carrier]
    if which == "capacity":
        df = stos.eval("p_nom_opt * max_hours").groupby(groups).sum()
    elif which == "operation":
        df = (
            (weights @ n.storage_units_t.state_of_charge[stos.index])
            .groupby(groups)
            .sum()
        )
    res.append(df)

    carriers = specs.get("Store", [])
    assert_carriers_existent(n, carriers, "Store")
    sus = n.stores.query("carrier in @carriers")
    groups = [sus.bus.map(location), sus.carrier]
    if which == "capacity":
        df = sus.groupby(groups).e_nom_opt.sum()
    elif which == "operation":
        df = (weights @ n.stores_t.e[sus.index]).groupby(groups).sum()
    res.append(df)

    return pd.concat(res)


def get_carrier_transport(n, kind, config, which="capacity"):
    """
    Get the transport of a certain kind of carrier.

    The corresponding values are calculated on the basis of the config file entry
    'kind_to_carrier'.

    Parameters
    ----------
    n : pypsa.Network
    kind : str
    config : dict
    which: str, optional
        Whether to return the capacity or the operation of the carriers, by default 'capacity'
    """
    weights = n.snapshot_weightings.generators
    specs = config["constants"]["kind_to_carrier"].get(kind, {}).get("transport", {})
    res = {}

    carriers = specs.get("Link", [])
    assert_carriers_existent(n, carriers, "Link")
    links = n.links.query("carrier in @carriers")
    if which == "capacity":
        df = links.p_nom_opt
    elif which == "operation":
        df = weights @ n.links_t.p0[links.index]
    elif which == "carrier":
        df = carriers
    res["Link"] = df

    carriers = specs.get("Line", [])
    assert_carriers_existent(n, carriers, "Line")
    lines = n.lines.query("carrier in @carriers")
    if which == "capacity":
        df = lines.s_nom_opt
    elif which == "operation":
        df = weights @ n.lines_t.p0[lines.index]
    elif which == "carrier":
        df = carriers
    res["Line"] = df

    return res


def get_generation_carriers(n):
    carriers = []
    for c in n.one_port_components | n.branch_components:
        if n.df(c).empty:
            continue
        cars = n.df(c).carrier
        if c == "Link":
            cars = cars[~get_transmission_links(n)]
        carriers.append(set(cars))
    return list(set.union(*carriers))


def sanitize_locations(n):
    if "EU" not in n.buses.index:
        n.add("Bus", "EU", x=-5.5, y=46)
        n.buses.loc["EU", "location"] = "EU"
        n.buses.loc["co2 atmosphere", "location"] = "EU"
    n.buses["x"] = n.buses.location.map(n.buses.x)
    n.buses["y"] = n.buses.location.map(n.buses.y)


def fill_missing_carriers(n):
    for c in n.iterate_components(n.one_port_components | n.branch_components):
        new_carriers = set(c.df.carrier.unique()) - set(n.carriers.index)
        if new_carriers:
            n.madd("Carrier", list(new_carriers), nice_name=list(new_carriers))


def modify_carrier_names(n):
    n.add("Carrier", "other", nice_name="Other", color="lightgrey")
    n.add("Carrier", "offwind", nice_name="Offshore Wind", color="#6895dd")
    n.mremove("Carrier", ["offwind-ac", "offwind-dc"])
    n.generators.loc[
        n.generators.carrier.str.startswith("offwind"), "carrier"
    ] = "offwind"
    replace = [f"(?i){s} " for s in SITES]
    for c in n.iterate_components(
        n.one_port_components | n.branch_components | {"Bus"}
    ):
        c.df.carrier.replace(replace, "", regex=True, inplace=True)
    n.carriers = n.carriers.sort_values("nice_name")
    n.carriers.index = n.carriers.index.to_series().replace(replace, "", regex=True)
    n.carriers.nice_name = n.carriers.nice_name.replace(replace, "", regex=True)
    n.carriers.nice_name = n.carriers.nice_name.replace(
        "solid biomass", "biomass", regex=True
    )
    n.carriers = n.carriers[~n.carriers.index.duplicated()]


def update_colors(n):
    config = yaml.load(open(root / "config" / "config.plotting.yaml"), yaml.CFullLoader)
    colors = pd.Series(config["plotting"]["tech_colors"])
    n.carriers.color.update(colors)


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
    # replace abbreviations with capital letters
    config = yaml.load(open(root / "config" / "config.plotting.yaml"), yaml.CFullLoader)
    nice_names = pd.Series(config["plotting"]["nice_names"])
    n.carriers.nice_name.update(nice_names)

    replace = {
        "H2": "H$_2$",
        "Hydrogen": "H$_2$",
        "Co2": "CO$_2$",
        "Chp": "CHP",
        "Dac": "DAC",
        "Smr": "SMR",
        " Cc": " CC",
        "Ocgt": "OCGT",
        "Ac": "AC",
        "Dc": "DC",
        "Sabatier": "Methanation",
        "Land Transport Ev": "Electric Vehicles",
        "Agriculture Machinery Oil Emissions": "Agriculture Emission",
        "Oil Emissions": "Aviation and Petrochemical\nEmission",
        "Shipping Methanol Emissions": "Shipping Methanol Emission",
        "Process Emissions": "Process Emission",
        "Allam": "Allam Cycle",
        "Sequestered": "Sequestration",
        "Gas For Industry": "Gas-based Industry Process",
        "Biomass For Industry": "Biomass-based Industry Process",
    }
    nice_names = n.carriers.nice_name.str.title()
    n.carriers.nice_name = nice_names.replace(replace, regex=True)


def add_carrier_groups(n):
    config = yaml.load(open(root / "config" / "config.plotting.yaml"), yaml.CFullLoader)
    groups = pd.Series(config["plotting"]["technology_groups"])
    colors = pd.Series(config["plotting"]["technology_group_colors"])
    n.carriers["group"] = groups.reindex(n.carriers.index, fill_value="")
    n.carriers["group_color"] = n.carriers.group.map(colors).fillna("")


def override_component_attrs(directory):
    attrs = Dict({k: v.copy() for k, v in component_attrs.items()})

    for component, list_name in components.list_name.items():
        fn = directory / f"{list_name}.csv"
        if fn.exists():
            overrides = pd.read_csv(fn, index_col=0, na_values="n/a")
            attrs[component] = overrides.combine_first(attrs[component])

    return attrs


def mock_snakemake(rulename, **wildcards):
    """
    This function is expected to be executed from the 'scripts'-directory of '
    the snakemake project. It returns a snakemake.script.Snakemake object,
    based on the Snakefile.
    If a rule has wildcards, you have to specify them in **wildcards.
    Parameters
    ----------
    rulename: str
        name of the rule for which the snakemake object should be generated
    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are
        needed.
    """
    import snakemake as sm
    import os
    from pypsa.descriptors import Dict
    from snakemake.script import Snakemake
    from packaging.version import Version, parse

    script_dir = Path(__file__).parent.resolve()
    assert (
        Path.cwd().resolve() == script_dir
    ), f"mock_snakemake has to be run from the repository scripts directory {script_dir}"
    try:
        os.chdir(script_dir.parent.parent)
        for p in sm.SNAKEFILE_CHOICES:
            if os.path.exists(p):
                snakefile = p
                break
        kwargs = (
            dict(rerun_triggers=[]) if parse(sm.__version__) > Version("7.7.0") else {}
        )
        workflow = sm.Workflow(snakefile, overwrite_configfiles=[], **kwargs)
        workflow.include(snakefile)
        workflow.global_resources = {}
        rule = workflow.get_rule(rulename)
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = Dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i in range(len(io)):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

        return snakemake

    finally:
        os.chdir(script_dir)


def print_kind_carriers(n, carriers):
    """Print all components and their carriers"""
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
                efficiency = f"efficiency{bus[-1] if int(bus[-1]) > 1 else ''}"
                # print(df.query("carrier == @e")[efficiency].mean())
