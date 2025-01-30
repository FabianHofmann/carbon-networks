import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import math
from common import (
    import_network,
    mock_snakemake,
    sort_rows_by_relative_diff,
    get_ordered_handles_labels,
)


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_bar_compare_diff",
        ext="png",
        clusters=90,
        comparison="emission-reduction-0.1",
    )

sns.set_theme(**snakemake.params["theme"])
labels = snakemake.config["labels"]

subtract = False
df = {}
carriers = []
for path in snakemake.input.networks:
    n = import_network(path, remove_gas_store_capital_cost=True)
    capex = n.statistics.capex()
    opex = n.statistics.opex(aggregate_time="sum")

    costs = capex.add(opex, fill_value=0)
    costs = costs.droplevel(0)[lambda x: x > 0]
    costs = costs.groupby(costs.index).sum()

    carriers.append(n.carriers)

    if not subtract:
        subtract = True
        ref = costs
        key = snakemake.params.labels[n.meta["wildcards"]["run"]]
    else:
        df[key] = ref.sub(costs, fill_value=0)
        subtract = False

diff = pd.concat(df, axis=1).fillna(0)
carriers = pd.concat(carriers).drop_duplicates()

norm = 1e9
unit = "bn€/a"

groups = snakemake.config["plotting"]["technology_groups"]
groups = {v: groups[k] for k, v in n.carriers.nice_name.items()}
# keep sequestration separate
groups["CO$_2$ Sequestration"] = "CO$_2$ Sequestration"
grouped = diff.mul(-1).groupby(groups).sum().div(norm)
# grouped = grouped[grouped.max(axis=1) > 5e4]
colors = snakemake.config["plotting"]["technology_group_colors"]
grouped = sort_rows_by_relative_diff(grouped)
rename = {"Carbon Capt. at Point Sources": "Carbon Capture\nat Point Sources"}
grouped = grouped.rename(rename, axis=0)
colors = {rename.get(k, k): v for k, v in colors.items()}

fig, ax = plt.subplots(
    1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
)

if not snakemake.config["configs"]["test"]:
    grouped = grouped[grouped.round(0).ne(0).any(axis=1)]


def rounded(x):
    # hard set this to figures in the appendix
    if abs(float(x)) <= 0.1:
        return ""
    elif grouped.abs().sum().sum() < 50:
        return round(x, 1)
    else:
        return int(round(x, 0))


grouped.T.plot(kind="bar", stacked=True, ax=ax, color=colors, legend=True, alpha=0.9)
for container in ax.containers:
    if abs(container.datavalues).sum() > grouped.abs().sum().sum() / 20:
        ax.bar_label(
            container,
            fmt=lambda x: rounded(x),
            label_type="center",
            fontsize=7,
            color="grey",
        )

pad = (abs(ax.get_ylim()[1]) + abs(ax.get_ylim()[0])) * 0.02
bbox = {"boxstyle": "circle", "facecolor": "none", "pad": 0.2, "edgecolor": "k"}
for i in range(len(grouped.columns)):
    col = grouped.iloc[:, i]
    if (val := col.sum()) != 0:
        val = rounded(val)
        y = col[col.ge(0)].sum() + pad
        ax.text(i, y, f"{val}", ha="center", va="bottom", fontsize=7, bbox=bbox)

ax.axhline(0, color="black", lw=0.5)
ax.set_ylabel(f"Net Investment Change [{unit}]")
ax.set_xlabel("")
# ax.grid(axis="y", alpha=0.5)

ax.legend().remove()
ax.legend(
    *get_ordered_handles_labels(ax, grouped),
    loc="center left",
    bbox_to_anchor=(1, 0.45),
    frameon=False,
    ncol=1,
)

sns.despine()
fig.savefig(snakemake.output.figure, dpi=300)
grouped.round(3).to_csv(snakemake.output.table)
