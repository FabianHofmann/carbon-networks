import os
import cartopy
from cartopy import crs as ccrs
from common import (
    import_network,
    mock_snakemake,
)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import geopandas as gpd


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_sequestration_map",
        run="half-price",
        clusters=40,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])

crs = ccrs.EqualEarth()
n = import_network(snakemake.input.network)
offshore_regions = gpd.read_file(
    snakemake.input.offshore_regions, index_col=0
).set_index("name")
offshore_regions = offshore_regions.to_crs(crs.proj4_init)
offshore_regions["potential"] = (
    n.stores.set_index("location").query("carrier == 'co2 stored'").e_nom_max.div(1e6)
)  # Mt/a


fig, ax = plt.subplots(
    figsize=snakemake.params.settings["figsize"], subplot_kw={"projection": crs}
)
offshore_regions.plot(
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
# draw costal lines
ax.coastlines(resolution="10m", color="black", linewidth=0.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
ax.set_aspect("equal")
ax.axis("off")
ax.set_extent(offshore_regions.total_bounds[[0, 2, 1, 3]], crs=crs)

fig.savefig(snakemake.output.map, bbox_inches="tight")
