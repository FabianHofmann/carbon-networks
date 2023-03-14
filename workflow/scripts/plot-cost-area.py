import os
from pathlib import Path
import warnings
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from common import (
    import_network,
    mock_snakemake,
)

sns.set_theme(style="white", context="paper", rc={"patch.linewidth": 0.1}, font="serif")

alpha = 1
region_alpha = 0.8


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_cost_area",
        design="co2network",
    )


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    capex = n.statistics.capex()
    opex = n.statistics.opex(aggregate_time="sum")

    costs = capex.add(opex, fill_value=0)
    costs = costs.droplevel(0)

    design, sequestration = Path(path).stem.split("_")
    key = (snakemake.config["labels"][design], int(sequestration))

    df[key] = costs
df = pd.concat(df, axis=1)


# do some grouping


norm = 1e6
unit = "M€"
nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)
sort_by_color = (
    lambda df: df.assign(color=colors[df.index])
    .sort_values(by="color")
    .drop("color", axis=1)
)

fig, ax = plt.subplots(1, 1, figsize=(3, 3.5), layout="constrained")
df.plot(kind="area", stacked=True, ax=ax, color=colors.to_dict(), alpha=0.8, lw=0)

ax.axhline(0, color="k", lw=1)
ax.set_xlabel("Sequestration Potential [Mt]")
ax.set_ylabel("System cost [{unit}]")
ax.grid(axis="y", alpha=0.5)

sns.despine()

fig.savefig(snakemake.output[0], bbox_inches="tight", dpi=300)
