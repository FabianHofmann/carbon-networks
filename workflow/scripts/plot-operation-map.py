import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from common import import_network, mock_snakemake

from plotting import (
    create_carrier_network,
    plot_map,
    add_legend,
    get_carrier_network_plotting_data,
)
import geopandas as gpd


sns.set_context("paper")
column = "Optimal Capacity"
alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():

    snakemake = mock_snakemake(
        "plot_operation_map",
        simpl="",
        lv=1.2,
        clusters=181,
        opts="",
        sector_opts="Co2L0-365H-T-H-B-I-A-solar+p3-linemaxext15-seq200",
        planning_horizons=2050,
        kind="electricity",
    )

n = import_network(snakemake.input.network)
offshore_regions = gpd.read_file(snakemake.input.offshore_regions).set_index("name")
onshore_regions = gpd.read_file(snakemake.input.onshore_regions).set_index("name")
kind = snakemake.wildcards.kind


fig, axes = plt.subplots(
    1, 2, figsize=(9, 4), subplot_kw={"projection": ccrs.EqualEarth()}
)

carriers = snakemake.config["constants"]["carrier_to_buses"][kind]
o = create_carrier_network(n, kind, carriers, include_eu=True)
data = get_carrier_network_plotting_data(o, "Operation")
specs = snakemake.config["plotting"]["operation_map"][kind]
bus_sizes = data.bus_sizes.copy()

if kind == "carbon":
    regions = offshore_regions
else:
    regions = onshore_regions

ax = axes[0]
data.bus_sizes = bus_sizes[bus_sizes > 0]
plot_map(
    ax,
    o,
    regions,
    data,
    bus_scale=float(specs["bus_scale"]),
    branch_scale=float(specs["branch_scale"]),
    alpha=alpha,
    branch_alpha=alpha,
    region_alpha=region_alpha,
    region_cmap=specs["region_cmap"],
    region_unit=specs["region_unit"],
)
fig.canvas.draw()
add_legend(
    ax,
    o,
    bus_scale=float(specs["bus_scale"]),
    branch_scale=float(specs["branch_scale"]),
    bus_sizes=specs["bus_sizes"],
    branch_sizes=specs["branch_sizes"],
    alpha=alpha,
    gen_carriers=o.carriers.loc[data.bus_sizes.index.unique(1)],
)
ax.set_extent(regions.total_bounds[[0, 2, 1, 3]])
ax.set_title(kind.title() + " Production")

data = get_carrier_network_plotting_data(o, "Operation")
data.bus_sizes = -bus_sizes[bus_sizes < 0]
ax = axes[1]
plot_map(
    ax,
    o,
    regions,
    data,
    bus_scale=float(specs["bus_scale"]),
    branch_scale=float(specs["branch_scale"]),
    alpha=alpha,
    branch_alpha=alpha,
    region_alpha=region_alpha,
    region_cmap=specs["region_cmap"],
    region_unit=specs["region_unit"],
)
fig.canvas.draw()
add_legend(
    ax,
    o,
    bus_scale=float(specs["bus_scale"]),
    branch_scale=float(specs["branch_scale"]),
    bus_sizes=specs["bus_sizes"],
    branch_sizes=specs["branch_sizes"],
    alpha=alpha,
    gen_carriers=o.carriers.loc[data.bus_sizes.index.unique(1)],
)
ax.set_extent(regions.total_bounds[[0, 2, 1, 3]])
ax.set_title(kind.title() + " Consumption")

fig.tight_layout()
fig.savefig(snakemake.output.map, bbox_inches="tight", dpi=300, facecolor="white")
