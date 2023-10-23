import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pypsa.plot import add_legend_circles, add_legend_patches, add_legend_lines
from common import (
    import_network,
    mock_snakemake,
    get_carrier_storage,
    get_carrier_transport,
    get_carrier_production,
)
import geopandas as gpd


column = "Optimal Capacity"
alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_capacity_country_bar",
        run="baseline",
        clusters=90,
        ext="pdf",
    )

sns.set_theme(**snakemake.params["theme"])

n = import_network(snakemake.input.network)
n.buses.country = n.buses.country.where(
    n.buses.country.ne(""), n.buses.location.map(n.buses.country)
)
config = snakemake.config
which = "capacity"
colors = n.carriers.set_index("nice_name").color.to_dict()

s = n.statistics

ds = s.optimal_capacity(groupby=s.get_country_and_carrier)
df = ds.Generator.unstack(fill_value=0).filter(regex="(?i)wind|solar|biomass")
df = df / 1e3
# df = ds.groupby(level=[1,2]).sum().unstack(fill_value=0).T

# %%
df.plot.bar(color=colors, subplots=True, legend=False, figsize=(5, 7), grid=True)
plt.tight_layout()
sns.despine()
plt.xlabel("")
plt.ylabel("GW")

plt.savefig(snakemake.output[0], bbox_inches="tight")
# %%
