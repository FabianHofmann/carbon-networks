import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from common import (
    get_ordered_handles_labels,
    import_network,
    mock_snakemake,
    sort_rows_by_relative_diff,
)

alpha = 1
region_alpha = 0.8

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_bar", ext="png", clusters=90, comparison="sequestration"
    )

sns.set_theme(**snakemake.params["theme"])

df = {}
objective = {}
for path in snakemake.input.networks:
    n = import_network(path, remove_gas_store_capital_cost=True, revert_dac=False)
    capex = n.statistics.capex()
    opex = n.statistics.opex(aggregate_time="sum")

    costs = capex.add(opex, fill_value=0)
    costs = costs.droplevel(0)[lambda x: x > 0]
    costs = costs.groupby(costs.index).sum()

    run = n.meta["wildcards"]["run"]
    labels = snakemake.params.labels
    key = labels[run]
    if (comparison := snakemake.wildcards["comparison"]) in labels["overwrites"]:
        key = labels["overwrites"][comparison].get(run, key)

    df[key] = costs
    objective[key] = n.objective

    # with gas storage capex removed, the following should hold
    # total_cost = costs.sum() - n.statistics.installed_capex().sum()
    # assert total_cost.round(0) / 1e9 == n.objective.round(0) / 1e9

df = pd.concat(df, axis=1)
objective = pd.Series(objective)

groups = snakemake.config["plotting"]["technology_groups"]
groups = {v: groups[k] for k, v in n.carriers.nice_name.items()}

grouped = df.groupby(groups).sum()
grouped = grouped[grouped.max(axis=1) > 5e4]

colors = snakemake.config["plotting"]["technology_group_colors"]

norm = 1e9
unit = "bn€/a"
grouped = sort_rows_by_relative_diff(grouped).div(norm)
rename = {"Carbon Capt. at Point Sources": "Carbon Capture\nat Point Sources"}
grouped = grouped.rename(rename, axis=0)
colors = {rename.get(k, k): v for k, v in colors.items()}
settings = snakemake.params.settings
overwrites = settings.pop("overwrites", {}).get(snakemake.wildcards.comparison, {})
settings.update(overwrites)

fig, ax = plt.subplots(1, 1, figsize=settings["figsize"], layout="constrained")

defaults = dict(kind="bar", stacked=True, rot=90, lw=0.2, alpha=0.8)
kwargs = {**defaults, **settings.get("kwargs", {})}
grouped.T.plot(ax=ax, color=colors, **kwargs)
for container in ax.containers:
    if container.datavalues.sum() > grouped.sum().sum() / 20:
        ax.bar_label(
            container,
            label_type="center",
            fmt=lambda x: int(round(x, 0)),
            fontsize=7,
            color="grey",
        )

pad = 5
total = grouped.sum()
baseline = total[total.index[0]]
decrease = 1 - total / baseline
for i, (val, y) in enumerate(zip(decrease, total)):
    if val == 0:
        ax.text(
            i,
            y + pad,
            f"{baseline:.0f}",
            ha="center",
            va="bottom",
            fontsize=7,
            color="gray",
        )
    else:
        ax.text(
            i,
            y + pad,
            f"{-val:+.1%}",
            ha="center",
            va="bottom",
            fontsize=7,
            color="gray",
        )


ax.axhline(0, color="k", lw=1)
ax.set_ylabel(f"System cost [{unit}]")

handles, labels = get_ordered_handles_labels(ax, grouped, wrap=22)
ax.legend(
    handles,
    labels,
    loc="center left",
    bbox_to_anchor=(1, 0.45),
    frameon=False,
    ncol=1,
)

sns.despine()
fig.savefig(snakemake.output[0], dpi=300)
grouped.fillna(0).round(3).to_csv(snakemake.output.table)
