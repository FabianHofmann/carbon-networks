import os
from pathlib import Path
from common import (
    import_network,
    mock_snakemake,
)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_captureshare_line",
        clusters=40,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])

config = snakemake.config
which = "capacity"


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    capacity = n.statistics.optimal_capacity()
    capacity = capacity.droplevel(0)[lambda x: x > 0]

    switch_techs = capacity.index[capacity.index.str.contains(" CC")]
    caps = pd.DataFrame(
        {
            "cc": capacity[switch_techs].rename(lambda ds: ds[:-3]),
            "original": capacity[switch_techs.str[:-3]],
        }
    )
    caps = caps["cc"] / caps.sum(1)

    key = snakemake.params.labels[n.meta["wildcards"]["run"]]
    df[key] = caps

df = pd.concat(df, axis=1)

nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)

fig, ax = plt.subplots(
    figsize=snakemake.params.settings["figsize"],
    # layout="constrained",
    sharex=True,
)

df.mul(100).T.plot(
    kind="line",
    ax=ax,
    color=colors.to_dict(),
    alpha=0.8,
    rot=90,
)
ax.set_ylabel(f"CC share [%]")
ax.grid(axis="y", alpha=0.5)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)

sns.despine()

fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
