import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from common import (
    import_network,
    mock_snakemake,
    sort_rows_by_diff,
    get_ordered_handles_labels,
)

alpha = 1
region_alpha = 0.8

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_compare_diff_bar", ext="png", clusters=90, comparison="baseline"
    )

sns.set_theme(**snakemake.params["theme"])
labels = snakemake.config["labels"]

subtract = False
df = {}
carriers = []
for path in snakemake.input.networks:
    n = import_network(path)
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
groups = {v: groups[k] for k, v in n.carriers.nice_name.drop("").items()}
grouped = diff.mul(-1).groupby(groups).sum().div(norm)
# grouped = grouped[grouped.max(axis=1) > 5e4]
colors = snakemake.config["plotting"]["technology_group_colors"]
grouped = sort_rows_by_diff(grouped)
rename = {"Carbon Capt. at Point Sources": "Carbon Capture\nat Point Sources"}
grouped = grouped.rename(rename, axis=0)
colors = {rename.get(k, k): v for k, v in colors.items()}

fig, ax = plt.subplots(
    1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
)

grouped = grouped[grouped.round(0).ne(0).any(axis=1)]
grouped.T.plot(
    kind="bar", stacked=True, ax=ax, color=colors, lw=0, legend=True, alpha=0.8
)

# ax.vlines(0, -0.5, len(diff) - 0.5, color="k", lw=0.5)

ax.set_ylabel(f"Net Investment Change [{unit}]")
ax.set_xlabel("")
ax.grid(axis="y", alpha=0.5)

ax.legend().remove()
ax.legend(
    *get_ordered_handles_labels(ax, grouped),
    loc="center left",
    bbox_to_anchor=(1, 0.45),
    frameon=False,
    ncol=1,
)

sns.despine()
fig.savefig(snakemake.output[0], dpi=300)
