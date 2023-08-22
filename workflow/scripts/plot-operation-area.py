import os
from pathlib import Path
from common import (
    import_network,
    mock_snakemake,
    get_carrier_consumption,
    get_carrier_production,
    group_small_contributions,
    sort_rows_by_diff,
)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


sns.set_theme(**snakemake.params["theme"])


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_operation_area", sigma="co2network", kind="carbon", ext="png"
    )


kind = snakemake.wildcards.kind
config = snakemake.config
labels = config["labels"]
which = "operation"


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    prod = get_carrier_production(n, kind, config, which).groupby(level=1).sum()
    cons = get_carrier_consumption(n, kind, config, which).groupby(level=1).sum()

    if (prod.sum() / cons.sum()).round(3) != 1:
        print(f"Warning: {kind} production and consumption are not equal.")

    design, sequestration = Path(path).stem.split("_")
    # key = (labels[design], int(sequestration))
    key = int(sequestration)

    df[key] = pd.concat([prod, cons], keys=["Production", "Consumption"])
df = pd.concat(df, axis=1)


norm = 1e6
unit = "TWh" if kind not in ["carbon", "co2"] else "Mt"
nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)

production = (
    df.loc["Production"].rename(index=nice_name).div(norm).pipe(sort_rows_by_diff).T
)
consumption = (
    df.loc["Consumption"].rename(index=nice_name).div(norm).pipe(sort_rows_by_diff).T
)

production = group_small_contributions(production, 1)
consumption = group_small_contributions(consumption, 1)

fig, ax = plt.subplots(
    1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
)

production.plot(
    kind="area",
    stacked=True,
    ax=ax,
    color=colors.to_dict(),
    alpha=0.8,
    lw=0,
)
h, l = ax.get_legend_handles_labels()
ax.legend().remove()
legend = fig.legend(
    h[::-1],
    l[::-1],
    loc="upper left",
    bbox_to_anchor=(1, 0.95),
    frameon=False,
    ncol=1,
    title="Production",
    labelspacing=0.3,
)
fig.add_artist(legend)

consumption.mul(-1).plot(
    kind="area", stacked=True, ax=ax, color=colors.to_dict(), alpha=0.8, lw=0
)
h, l = ax.get_legend_handles_labels()
handle_map = (
    pd.Series(h, index=l)[lambda ds: ~ds.index.duplicated()]
    .reindex(consumption.columns)
    .dropna()
)
h = handle_map.values
l = handle_map.index

ax.legend().remove()
legend = fig.legend(
    h,
    l,
    loc="lower left",
    bbox_to_anchor=(1, 0.05),
    frameon=False,
    ncol=1,
    title="Consumption",
    labelspacing=0.3,
)
fig.add_artist(legend)

ax.axhline(0, color="k", lw=1)
ax.set_xlim(production.index.min(), production.index.max())
ax.set_ylim(-consumption.sum(1).max() * 1.1, production.sum(1).max() * 1.1)
ax.set_xlabel("Sequestration Potential [Mt]")
ax.set_ylabel(f"{labels[kind]} [{unit}]")
if snakemake.params.settings.get("title", True):
    ax.set_title(f"{labels[kind]} Balance {labels[design]}")
ax.grid(axis="y", alpha=0.5)

sns.despine()

# fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
