import os
import re
from pathlib import Path
from common import (
    import_network,
    mock_snakemake,
)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_captureshare_line",
        comparison="default",
        clusters=90,
        ext="png",
    )

sns.set_theme(**snakemake.params["theme"])

config = snakemake.config
which = "capacity"


share = {}
cf = {}
for path in snakemake.input.networks:
    n = import_network(path)
    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    s = n.statistics
    capacity = s.optimal_capacity()
    capacity = capacity.Link[lambda x: x > 0]

    switch_techs = capacity.index[capacity.index.str.contains(" CC")]
    caps = pd.DataFrame(
        {
            "cc": capacity[switch_techs].rename(lambda ds: ds[:-3]),
            "original": capacity[switch_techs.str[:-3]],
        }
    )
    caps = caps["cc"] / caps.sum(1)
    share[key] = caps.rename_axis("Carrier")

    def groupby(n, c: str, nice_names) -> pd.Series:
        return (
            n.df(c).carrier.replace(" CC", "", regex=True).replace(n.carriers.nice_name)
        )

    cf[key] = s.capacity_factor(groupby=groupby).Link[caps.index].rename_axis("Carrier")

share = pd.concat(share, axis=1, names="Model")
cf = pd.concat(cf, axis=1, names="Model")

data = pd.concat(
    [df.mul(100).stack() for df in [cf, share]],
    axis=1,
    keys=["Capacity Factor", "CC Share [%]"],
).reset_index()
# %%
fig, ax = plt.subplots(
    figsize=snakemake.params.settings["figsize"],
    layout="constrained",
)

nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)


plot = sns.scatterplot(
    ax=ax,
    data=data,
    x="Model",
    y="CC Share [%]",
    size="Capacity Factor",
    hue="Carrier",
    # style="Carrier",
    palette=colors.to_dict(),
    legend="auto",
)
handles, labels = ax.get_legend_handles_labels()
labels = [l + " %" if re.match(r"^\d+$", l) else l for l in labels]

pos = labels.index("Capacity Factor")
handles.insert(pos, Patch(facecolor="none", edgecolor="none", alpha=0))
labels.insert(pos, "")

# Add the invisible handles and labels to the legend
ax.legend(
    handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, ncol=1
)
ax.grid(axis="y", alpha=0.5)
ax.set_xlabel("")

plt.xticks(rotation=90)

sns.despine()

fig.savefig(snakemake.output[0], dpi=300)
