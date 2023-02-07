KINDS = [
    "electricity",
    "gas",
    "hydrogen",
    "carbon",
]
ONEPORT_COMPS = ["Generator", "Store", "StorageUnit"]
NODE_COMPS = ["Generator", "Link", "StorageUnit"]
REGION_COMPS = ["Store"]
BRANCH_COMPS = ["Line", "Link"]
COMPS = set(NODE_COMPS + REGION_COMPS + BRANCH_COMPS)
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

PLOT_SPECS = {
    "carbon": {
        "region_cmap": "Reds",
        "bus_scale": 1e-4,
        "branch_scale": 5e-3,
        "bus_sizes": [5_000, 10_000],
        "branch_sizes": [1_000, 2_000],
    },
    "hydrogen": {
        "region_cmap": "Blues",
        "bus_scale": 4e-5,
        "branch_scale": 3e-3,
        "bus_sizes": [5_000, 10_000],
        "branch_sizes": [1_000, 2_000],
    },
    "electricity": {
        "region_cmap": "Greens",
        "bus_scale": 5e-6,
        "branch_scale": 2e-4,
        "bus_sizes": [20_000, 50_000],
        "branch_sizes": [5_000, 10_000],
    },
    "gas": {
        "region_cmap": "Purples",
        "bus_scale": 2e-6,
        "branch_scale": 1e-4,
        "bus_sizes": [100_000, 200_000],
        "branch_sizes": [10_000, 20_000],
    },
}
