import os
from pathlib import Path
import cartopy.crs as ccrs
import seaborn as sns
import matplotlib.pyplot as plt
from pypsa.plot import add_legend_circles, add_legend_patches, add_legend_lines
import pandas as pd
from common import (
    import_network,
    mock_snakemake,
)

alpha = 1
region_alpha = 0.8

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_diff_map", ext="pdf", clusters=90, difference="default"
    )

sns.set_theme(**snakemake.params["theme"])
config = snakemake.config
labels = config["labels"]
specs = config["plotting"]["capacity_diff_map"]

df = {}
carriers = []
for path in snakemake.input.networks:
    n = import_network(path)
    s = n.statistics
    costs = s.capex(groupby=s.get_bus_and_carrier) + s.opex(
        groupby=s.get_bus_and_carrier
    )
    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    df[key] = costs
    carriers.append(n.carriers)


df = pd.concat(df, axis=1).fillna(0)
carriers = pd.concat(carriers).drop_duplicates()
colors = carriers.set_index("nice_name").color

norm = 1e9
unit = "bn€/a"

diff = (df[df.columns[0]] - df[df.columns[1]]).div(norm)
diff = diff.rename(n.buses.location, level=1)
diff = diff[diff.round(0) != 0].droplevel(0)
drop = ["H$_2$ Pipeline", "CO$_2$ Pipeline", "AC", "Heat Waste", "Gas Pipeline"]
diff = diff.drop(drop, level=1, errors="ignore")

pos = diff[diff > 0]
neg = diff[diff < 0]


fig, axes = plt.subplots(
    1,
    2,
    figsize=snakemake.params.settings["figsize"],
    squeeze=False,
    subplot_kw={"projection": ccrs.EqualEarth()},
)

bus_scale = float(specs["bus_scale"])
alpha = 0.8

for ax, ds, col in zip(axes.flat, [pos, neg], df.columns):
    n.plot(
        bus_sizes=ds.abs() * bus_scale,
        bus_alpha=alpha,
        ax=ax,
        bus_colors=colors,
        link_widths=0,
        line_widths=0,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
        margin=0.2,
    )
    ax.set_title(f"Higher Spendings {col}")

    legend_bus_sizes = specs["bus_sizes"]
    if legend_bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_scale for s in legend_bus_sizes],
            [f"{s} {unit}" for s in legend_bus_sizes],
            legend_kw={"bbox_to_anchor": (0, 1), "loc": "upper left", "frameon": False},
        )

gen_carriers = (
    carriers.set_index("nice_name").loc[diff.index.unique(1)].sort_values("color")
)
add_legend_patches(
    fig,
    gen_carriers.color,
    gen_carriers.index,
    patch_kw={"alpha": alpha},
    legend_kw={
        "bbox_to_anchor": (0.5, 0),
        "ncol": 5,
        "loc": "upper center",
        "frameon": False,
    },
)

ax.set_extent(snakemake.config["plotting"]["extent"])

# fig.tight_layout()
fig.savefig(
    snakemake.output[0],
    dpi=300,
    bbox_inches="tight",
)
