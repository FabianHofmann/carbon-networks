import os
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from common import (
    import_network,
    mock_snakemake,
    sort_rows_by_diff,
)

sns.set_theme(style="white", context="paper", rc={"patch.linewidth": 0.1}, font="serif")

alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake("plot_cost_area", design="co2network", ext="png")


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    capex = n.statistics.capex()
    opex = n.statistics.opex(aggregate_time="sum")

    costs = capex.add(opex, fill_value=0)
    costs = costs.droplevel(0)[lambda x: x > 0]
    costs = costs.groupby(costs.index).sum()

    design, sequestration = Path(path).stem.split("_")
    key = design, int(sequestration)

    df[key] = costs
df = pd.concat(df, axis=1)

groups = snakemake.config["plotting"]["technology_groups"]
grouped = df.groupby(groups).sum()
grouped = grouped[grouped.max(axis=1) > 5e4]

colors = snakemake.config["plotting"]["technology_group_colors"]

norm = 1e9
unit = "bn€/a"
sort_by_color = (
    lambda df: df.assign(color=colors).sort_values(by="color").drop("color", axis=1)
)
grouped = sort_rows_by_diff(grouped).div(norm)

fig, axes = plt.subplots(1, 2, figsize=(7, 3.0), layout="constrained")

for key, ax in zip(snakemake.config["scenarios"]["design"], axes):

    grouped[key].T.plot(kind="area", stacked=True, ax=ax, color=colors, alpha=0.8, lw=0)

    ax.axhline(0, color="k", lw=1)
    ax.set_xlabel("Sequestration Potential [Mt]")
    ax.set_ylabel(f"System cost [{unit}]")
    ax.grid(axis="y", alpha=0.5)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend().remove()
    ax.set_title(snakemake.config["labels"][key].title().replace("Co", "CO"))


sns.despine()
fig.legend(
    handles[::-1],
    labels[::-1],
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
    ncol=1,
)
fig.savefig(snakemake.output[0], bbox_inches="tight", dpi=300)
