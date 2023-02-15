#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 13:55:01 2023

@author: fabian
"""
import os
from pathlib import Path
import warnings
import pypsa
from pypsa.descriptors import Dict
from pypsa.components import components, component_attrs
import pandas as pd
import yaml
import numpy as np

root = Path(__file__).parent.parent.parent.resolve()
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


def get_transmission_links(n):
    # only choose transmission links
    return (
        (n.links.bus0.map(n.buses.location) != n.links.bus1.map(n.buses.location))
        & ~n.links.bus0.map(n.buses.location).str.contains("EU")
        & ~n.links.bus1.map(n.buses.location).str.contains("EU")
    )


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


def assign_location(n):
    for c in n.one_port_components | n.branch_components:
        df = n.df(c)

        if "location" not in df:
            df["location"] = np.nan

        ifind = pd.Series(df.index.str.find(" ", start=4), df.index)
        for i in ifind.value_counts().index:
            # these have already been assigned defaults
            # if i == -1:
            #     continue
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


def add_colors(n):
    config = yaml.load(
        open(root / "config" / "config.pypsa-eur-sec.yaml"), yaml.CFullLoader
    )
    colors = config["plotting"]["tech_colors"]
    n.carriers.color = n.carriers.color.where(n.carriers.color != "")
    n.carriers.color.update(colors)
    mapper = get_carrier_mapper(n)
    n.carriers.color = n.carriers.color.fillna(mapper.map(colors))


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
