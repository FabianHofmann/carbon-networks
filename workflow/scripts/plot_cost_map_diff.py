import os
import cartopy.crs as ccrs
import geopandas as gpd
import seaborn as sns
import matplotlib.pyplot as plt
from pypsa.plot import add_legend_circles, add_legend_patches, add_legend_lines
from pypsa.statistics import get_transmission_branches, get_transmission_carriers
import pandas as pd
from common import (
    import_network,
    mock_snakemake,
)

alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_map_diff",
        ext="png",
        clusters=90,
        difference="emission-reduction-0.1-full",
    )

sns.set_theme(**snakemake.params["theme"])
plt.rc("patch", linewidth=0.1)

config = snakemake.config
labels = config["labels"]
specs = config["plotting"]["cost_map_diff"]
cutoff_bus = specs["cutoff_bus"]
cutoff_branch = specs["cutoff_branch"]

networks = [
    import_network(path, remove_gas_store_capital_cost=True)
    for path in snakemake.input.networks
]
regions = gpd.read_file(snakemake.input.regions).set_index("name")
df_bus = {}
df_branch = {}
carriers = []
for n in networks:
    s = n.statistics

    # TODO: ensure this weird adjustment is not necessary anymore
    is_dac = n.links.carrier == "DAC"
    n.links.loc[is_dac, ["bus0", "bus1"]] = n.links.loc[is_dac, ["bus1", "bus0"]].values

    capex = s.capex(groupby=s.get_bus_and_carrier)
    opex = s.opex(aggregate_time="sum", groupby=s.get_bus_and_carrier)
    costs = capex.add(opex, fill_value=0)

    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    df_bus[key] = costs

    branches = get_transmission_branches(n)
    df_branch[key] = s.capex(groupby=False).loc[branches]

    carriers.append(n.carriers)


df_bus = pd.concat(df_bus, axis=1).fillna(0)
df_branch = pd.concat(df_branch, axis=1).fillna(0)
carriers = pd.concat(carriers).drop_duplicates()
colors = carriers.set_index("nice_name").color

norm_bus = 1e9
unit_bus = "bn€/a"
norm_branch = 1e6
unit_branch = "M€/a"

diff_bus = df_bus.diff(axis=1)[df_bus.columns[1]]
diff_bus = diff_bus[diff_bus.abs() > cutoff_bus]
diff_bus = diff_bus.div(norm_bus)
diff_bus = diff_bus.rename(n.buses.location, level=1)
diff_bus = diff_bus.droplevel(0)
drop = ["H$_2$ Pipeline", "CO$_2$ Pipeline", "AC", "Heat Waste", "Gas Pipeline"]
diff_bus = diff_bus.drop(drop, level="carrier", errors="ignore")

pos_bus = diff_bus[diff_bus > 0].abs()
neg_bus = diff_bus[diff_bus < 0].abs()

diff_branch = df_branch.diff(axis=1)[df_branch.columns[1]]
diff_branch = diff_branch[diff_branch.abs() > cutoff_branch].div(norm_branch)
pos_branch = diff_branch[diff_branch > 0].abs()
neg_branch = diff_branch[diff_branch < 0].abs()

fig, axes = plt.subplots(
    1,
    2,
    figsize=snakemake.params.settings["figsize"],
    squeeze=False,
    subplot_kw={"projection": ccrs.EqualEarth()},
)
bus_scale = float(specs["bus_scale"])
branch_scale = float(specs["branch_scale"])
alpha = 1

for ax, ds_bus, ds_branch, col, n in zip(
    axes.flat, [neg_bus, pos_bus], [neg_branch, pos_branch], df_bus.columns, networks
):
    links = n.links.index.intersection(ds_branch.Link.index)
    link_widths = ds_branch.Link[links]
    link_colors = n.links.carrier[links].map(n.carriers.color)

    if "Line" in ds_branch.index.unique(0):
        lines = n.lines.index.intersection(ds_branch.Line.index)
        line_widths = ds_branch.Line[lines]
        line_colors = n.lines.carrier[lines].map(n.carriers.color)
    else:
        line_widths = pd.Series()
        line_colors = pd.Series()

    n.plot(
        bus_sizes=ds_bus * bus_scale,
        bus_alpha=alpha,
        ax=ax,
        bus_colors=colors,
        link_widths=link_widths * branch_scale,
        link_colors=link_colors,
        link_alpha=alpha,
        line_widths=line_widths * branch_scale,
        line_colors=line_colors,
        line_alpha=alpha,
        color_geomap={"border": "darkgrey", "coastline": "darkgrey"},
        geomap="10m",
        boundaries=snakemake.config["plotting"]["extent"],
    )
    regions.plot(
        ax=ax,
        facecolor="whitesmoke",
        edgecolor="darkgrey",
        linewidth=0,
        alpha=region_alpha,
        transform=ccrs.PlateCarree(),
        aspect="equal",
    )

    title = col.replace("\n", " ")
    ax.set_title(f"Higher Spendings {title}")

    legend_bus_sizes = specs["bus_sizes"]
    if legend_bus_sizes is not None:
        add_legend_circles(
            ax,
            [s * bus_scale for s in legend_bus_sizes],
            [f"{s} {unit_bus}" for s in legend_bus_sizes],
            legend_kw={
                "bbox_to_anchor": (0, 0.95),
                "loc": "upper left",
                "frameon": True,
            },
        )
    legend_branch_sizes = specs["branch_sizes"]
    if legend_branch_sizes is not None:
        add_legend_lines(
            ax,
            [s * branch_scale for s in legend_branch_sizes],
            [f"{s} {unit_branch}" for s in specs["branch_sizes"]],
            legend_kw={
                "bbox_to_anchor": (0, 0.82),
                "loc": "upper left",
                "frameon": True,
            },
        )


gen_carriers = (
    carriers.set_index("nice_name").loc[diff_bus.index.unique(1)].sort_values("color")
)
add_legend_patches(
    fig,
    gen_carriers.color,
    gen_carriers.index,
    patch_kw={"alpha": alpha},
    legend_kw={
        "bbox_to_anchor": (0.4, 0),
        "ncol": 4,
        "loc": "upper center",
        "frameon": True,
        "title": "Production, Storage & Conversion",
    },
)

branch_carriers = get_transmission_carriers(n).unique(1)
branch_carriers = carriers.loc[branch_carriers].sort_values("color")
add_legend_patches(
    fig,
    branch_carriers.color,
    branch_carriers.nice_name,
    patch_kw={"alpha": alpha},
    legend_kw={
        "bbox_to_anchor": (0.85, 0),
        "ncol": 1,
        "loc": "upper center",
        "frameon": True,
        "title": "Transmission",
    },
)


fig.tight_layout()
fig.savefig(
    snakemake.output[0],
    dpi=300,
    bbox_inches="tight",
)
