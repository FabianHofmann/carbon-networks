import os
from pathlib import Path
from common import (
    import_network,
    mock_snakemake,
)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


sns.set_theme(style="white", context="paper", font="serif")


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_captureshare_line",
        design="co2network",
        ext="png",
    )


config = snakemake.config
which = "capacity"


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    capacity = n.statistics.optimal_capacity()
    capacity = capacity.droplevel(0)[lambda x: x > 0]

    switch_techs = capacity.index[capacity.index.str.contains("\*")]
    caps = pd.DataFrame(
        {
            "cc": capacity[switch_techs].rename(lambda ds: ds[:-1]),
            "original": capacity[switch_techs.str[:-1]],
        }
    )
    caps = caps["cc"] / caps.sum(1)

    design, sequestration = Path(path).stem.split("_")
    key = (snakemake.config["labels"][design], int(sequestration))
    # key = int(sequestration)

    df[key] = caps
df = pd.concat(df, axis=1)

nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)

fig, axes = plt.subplots(2, 1, figsize=(5, 4), layout="constrained", sharex=True)

for ax, col in zip(axes, df.columns.unique(0)):
    df[col].mul(100).sort_values(by=200, ascending=False).T.plot(
        kind="line",
        ax=ax,
        color=colors.to_dict(),
        alpha=0.8,
    )
    ax.set_xlabel("Sequestration Potential [Mt]")
    ax.set_ylabel(f"CC share [%]")
    ax.set_title(col.title().replace("Co", "CO"))
    ax.grid(axis="y", alpha=0.5)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)

sns.despine()

fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
