import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from common import (
    import_network,
    mock_snakemake,
)

alpha = 1
region_alpha = 0.8

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_bar_diff", ext="png", clusters=90, difference="default"
    )

sns.set_theme(**snakemake.params["theme"])
labels = snakemake.config["labels"]
cutoff = snakemake.config["plotting"]["cost_bar_diff"]["cutoff"]

df = {}
carriers = []
for path in snakemake.input.networks:
    n = import_network(path)
    capex = n.statistics.capex()
    opex = n.statistics.opex(aggregate_time="sum")

    costs = capex.add(opex, fill_value=0)
    costs = costs.droplevel(0)[lambda x: x > 0]
    costs = costs.groupby(costs.index).sum()

    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    df[key] = costs
    carriers.append(n.carriers)

df = pd.concat(df, axis=1).fillna(0)
carriers = pd.concat(carriers).drop_duplicates()

norm = 1e9
unit = "bn€/a"

diff = df[df.columns[0]] - df[df.columns[1]]
diff = diff[diff.abs() > cutoff].div(norm)

fig, ax = plt.subplots(
    1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
)

diff = diff[diff.round(1) != 0].sort_values()
colors = carriers.set_index("nice_name").color[diff.index]
diff.plot(kind="barh", ax=ax, color=colors, lw=0, legend=False)
for container in ax.containers:
    if abs(container.datavalues).sum() > diff.abs().sum().sum() / 20:
        ax.bar_label(
            container,
            label_type="center",
            fontsize=7,
            color="grey",
            labels=diff.round(0).astype(int),
        )

ax.vlines(0, -0.5, len(diff) - 0.5, color="k", lw=0.5)
ax.set_xlabel(f"Cost benefit [{unit}]")
ax.set_ylabel("")
ax.grid(axis="both", alpha=0.5)
ax.set_title(f"Cost difference ({df.columns[0]} - {df.columns[1]})")

sns.despine(left=True)
fig.savefig(snakemake.output[0], dpi=300)
