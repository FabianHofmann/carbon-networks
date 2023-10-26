import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
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


grouper = [n.generators.bus.map(n.buses.country), n.generators.carrier]
p_nom_opt = n.generators.p_nom_opt.groupby(grouper).sum()
p_nom_max = n.generators.p_nom_max.groupby(grouper).sum()


filter = "(?i)wind|solar|biomass"
p_nom_opt = p_nom_opt.unstack(fill_value=0).filter(regex=filter).div(1e3)
p_nom_max = p_nom_max.unstack(fill_value=0).filter(regex=filter).div(1e3)

p_nom_opt = p_nom_opt.rename(columns=n.carriers.nice_name)
p_nom_max = p_nom_max.rename(columns=n.carriers.nice_name)


# %%

fig, axes = plt.subplots(len(p_nom_opt.columns), 1, figsize=(6, 10), sharex=False)

colors = sns.color_palette("bright", len(p_nom_opt.index))
for ax, col in zip(axes, p_nom_opt.columns):

    title = col
    if p_nom_max[col].replace(np.inf, np.nan).isna().all():
        title += " (no limit)"
    else:
        p_nom_max[col].plot.bar(legend=False, ax=ax, color=colors, alpha=0.35)

    p_nom_opt[col].plot.bar(legend=False, ax=ax, color=colors)

    ax.grid(axis="y")
    ax.set_ylabel("GW")
    ax.set_title(title)
    ax.set_xlabel("")

plt.tight_layout()
sns.despine()

fig.savefig(snakemake.output[0], bbox_inches="tight")

# %%
