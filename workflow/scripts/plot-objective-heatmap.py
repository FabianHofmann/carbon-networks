# %%

import os
import pandas as pd
import pypsa
import matplotlib.pyplot as plt
import seaborn as sns
from common import (
    import_network,
    mock_snakemake,
    sort_rows_by_diff,
)

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_objectice_heatmap",
        comparison="emission-reduction",
        ext="pdf",
        clusters=90,
    )

sns.set_theme(**snakemake.params["theme"])

config = snakemake.config


df = {}
for path in snakemake.input.networks:
    n = import_network(path)

    # key = snakemake.params.labels[n.meta["wildcards"]["run"]]
    key = n.meta["wildcards"]["run"]

    df[key] = n.statistics.capex().sum() + n.statistics.opex().sum()

df = pd.Series(df)

index = [l if len(l) == 2 else ["net-neutral"] + l for l in df.index.str.split("/")]
df.index = pd.MultiIndex.from_tuples(index)
df = df.unstack()

df = df.loc[["net-neutral", "net-negative-0.05", "net-negative-0.1"]][
    ["baseline", "co2-only", "h2-only", "full"]
]
df = df.div(1e9).round(0).astype(int)  # bn €

data = df.rename(config["labels"], axis=1).rename(config["labels"], axis=0).T

fig, ax = plt.subplots(1, 1, figsize=snakemake.params.settings["figsize"])
annot = data.applymap(lambda x: f"{x} bn €")
sns.heatmap(
    data,
    ax=ax,
    cmap="RdBu_r",
    cbar=False,
    annot=annot,
    fmt="",
    annot_kws={"fontsize": 9},
)
# fig.tight_layout()
fig.savefig(snakemake.output[0])
