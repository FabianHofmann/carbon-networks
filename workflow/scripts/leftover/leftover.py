CARRIER_TO_PLOT = {
    "carbon": {
        "node_carrier": [
            "process emissions CC",
            "solid biomass for industry CC",
            "SMR CC",
            "gas for industry CC",
            "solid biomass CHP CC",
            "gas CHP CC",
            "DAC",
        ],
        "region_carrier": ["co2 stored"],
        "branch_carrier": ["CO2 pipeline"],
    },
    "hydrogen": {
        "node_carrier": ["H2 Electrolysis"],
        "region_carrier": ["H2"],
        "branch_carrier": ["H2 pipeline"],
    },
    "electricity": {
        "node_carrier": [
            "solar",
            "offwind",
            "onwind",
            "ror",
            "hydro",
            "PHS",
            "OCGT",
            "H2 Fuel Cell",
            "CCGT",
            "geothermal",
            "biomass",
        ],
        "branch_carrier": ["AC", "DC"],
        "region_carrier": ["battery", "home battery", "Li ion"],
    },
    "gas": {
        "node_carrier": ["gas"],
        "branch_carrier": ["gas pipeline"],
        "region_carrier": ["gas", "biogas"],
    },
}


CARRIER_GROUPS = {
    "Battery Infrastructure": [
        "battery",
        "battery discharger",
        "battery charger",
        "Li ion",
        "home battery",
    ],
    "Fossil Carriers": ["coal", "CCGT", "lignite", "OCGT"],
    "Hydrogen Infrastructure": [
        "H2",
        "H2 fuel cell",
        "H2 electrolysis",
        "H2 Electrolysis",
    ],
    "Nuclear": ["nuclear"],
    "Offshore Wind": ["offwind"],
    "Onshore Wind": ["onwind"],
    "Other Renewables": [
        "geothermal",
        "biomass",
        "PHS",
        "hydro",
        "ror",
        "biogas",
        "solid biomass",
    ],
    "Solar": ["solar"],
    "Gas": ["gas"],
    "Transmission System": ["AC", "DC"],
    "Heat": [f"{site} heat" for site in SITES]
    + [
        "heat",
        "air heat pump",
        "ground heat pump",
        "biomass boiler",
        "gas boiler",
        "oil boiler",
        "wood boiler",
        "solar thermal",
    ],
    "Water Tanks": [f"{site} water tanks" for site in SITES],
    "Solar Thermal": [f"{site} solar thermal" for site in SITES],
    "Carbon Capture": [
        "solid biomass CHP CC",
        "gas CHP CC",
        "process emissions CC",
        "solid biomass for industry CC",
        "SMR CC",
        "gas for industry CC",
        "DAC",
    ],
    "Electricity": ["AC", "DC", "electricity distribution grid"],
}

CARRIER_GROUPS_r = {v: k for k, v in pd.Series(CARRIER_GROUPS).explode().items()}


def add_carrier_groups(n):
    n.carriers["group"] = get_carrier_mapper(n).map(CARRIER_GROUPS_r)


def get_carrier_network_data(n, kind):
    """
    Given a network, the function returns the following data:

    - name: Name of the data
    - bus_sizes: Series of bus sizes with MultiIndex (location, carrier)
    - branch_widths: Series of link widths
    - line_colors: Series of line colors
    - storage_capacities: Series of storage capacities per bus
    - carriers: Dataframe of carriers
    """
    specs = CARRIER_TO_PLOT[kind]
    kwargs = dict(aggregate_time="sum")

    groups = ["location", "carrier"]
    node = (
        n.statistics(NODE_COMPS, groupby=groups, **kwargs)
        .reindex(specs["node_carrier"], level="carrier")
        .groupby(level=groups)
        .sum()
    )

    region = (
        n.statistics(REGION_COMPS, groupby=groups, **kwargs)
        .reindex(specs["region_carrier"], level="carrier")
        .groupby(level="location")
        .sum()
    )

    branch_color = n.carriers.color[specs["branch_carrier"][0]]
    groups = ["interconnection", "carrier"]
    branch = (
        n.statistics(BRANCH_COMPS, groupby=groups, **kwargs)
        .reindex(specs["branch_carrier"], level="carrier")
        .groupby(level=groups)
        .sum()
    )

    # correct links by efficiency
    cols = ["Optimal Capacity"]
    efficiency = n.links.groupby("carrier").efficiency.mean()
    branch[cols] = branch[cols].div(efficiency, level="carrier", axis=0)
    branch = branch.groupby(level="interconnection").sum()

    carriers = specs["node_carrier"] + specs["region_carrier"] + specs["branch_carrier"]

    return Dict(
        {
            "node": node,
            "branch": branch,
            "branch_color": branch_color,
            "region": region,
            "carriers": carriers,
        }
    )


def create_dummy_network(n):
    dummy = pypsa.Network()
    buses = n.buses.loc[n.buses.location.unique(), ["x", "y"]]
    dummy.madd("Bus", buses.index, **buses)
    connections = (
        n.branches()
        .drop_duplicates(subset=["interconnection"])
        .dropna(subset=["interconnection"])
        .set_index("interconnection")
        .drop(columns=["bus0", "bus1"])
        .rename(columns={"location0": "bus0", "location1": "bus1"})[["bus0", "bus1"]]
    )
    dummy.madd("Link", connections.index, **connections)
    dummy.madd("Carrier", n.carriers.index, **n.carriers)
    return dummy


# rule create_gifs:
#     input:
#         # use input function here
#         expand(
#             "results/comparison-{comparison}/gifs/{clusters}_nodes/capacity_map_{kind}.gif",
#             kind=config["constants"]["kinds"],
#             clusters=config["scenario"]["clusters"],
#         ),
#         expand(
#             "results/comparison-{comparison}/gifs/{clusters}_nodes/balance_map_{kind}.gif",
#             kind=config["constants"]["kinds"],
#             clusters=config["scenario"]["clusters"],
#         ),
#     output:
#         "results/comparison-{comparison}/gifs/{clusters}_nodes/.gifs_created",
#     shell:
#         "touch {output}"
