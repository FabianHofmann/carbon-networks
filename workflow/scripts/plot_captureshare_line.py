import os
import re
import numpy as np
from pathlib import Path
from common import (
    import_network,
    mock_snakemake,
    groupby_carrier_across_cc,
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

    cf[key] = (
        s.capacity_factor(groupby=groupby_carrier_across_cc)
        .Link[caps.index]
        .rename_axis("Carrier")
    )

share = pd.concat(share, axis=1, names="Model")
cf = pd.concat(cf, axis=1, names="Model")

data = pd.concat(
    [df.mul(100).stack() for df in [cf, share]],
    axis=1,
    keys=["Capacity Factor", "Carbon Capture Share [%]"],
).reset_index()

fig, ax = plt.subplots(
    figsize=snakemake.params.settings["figsize"],
    layout="constrained",
)

nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)


model_to_number = {k: i for i, k in enumerate(data.Model.unique())}
data["model_jittered"] = (
    data["Model"].replace(model_to_number) + (np.random.rand(len(data)) - 0.5) * 0.3
)


plot = sns.scatterplot(
    ax=ax,
    data=data,
    x="model_jittered",
    y="Carbon Capture Share [%]",
    size="Capacity Factor",
    hue="Carrier",
    # style="Carrier",
    palette=colors.to_dict(),
    legend="auto",
    alpha=0.8,
    marker="o",
    sizes=(20, 100),
    linewidth=0.5,
)


handles, labels = ax.get_legend_handles_labels()
labels = [l + " %" if re.match(r"^\d+$", l) else l for l in labels]

pos = labels.index("Capacity Factor")
handles.insert(pos, Patch(facecolor="none", edgecolor="none", alpha=0))
labels.insert(pos, "")

# Add the invisible handles and labels to the legend
ax.legend(
    handles,
    labels,
    loc="center left",
    bbox_to_anchor=(1.0, 0.5),
    frameon=False,
    ncol=1,
)
ax.grid(axis="y", alpha=0.5)
ax.set_xlabel("")

plt.xticks(
    ticks=list(model_to_number.values()),
    labels=list(model_to_number.keys()),
    rotation=90,
)

sns.despine()

fig.savefig(snakemake.output.figure, dpi=300)
data.to_csv(snakemake.output.table)

# ALTERATIVE PLOT
# fig, ax = plt.subplots(
#     figsize=snakemake.params.settings["figsize"],
#     layout="constrained",
# )
# plot = sns.scatterplot(
#     ax=ax,
#     data=data,
#     style="Model",
#     y="CC Share [%]",
#     x="Capacity Factor",
#     hue="Carrier",
#     palette=colors.to_dict(),
#     legend="auto",
# )
# ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, ncol=1)
# x = data["Capacity Factor"]
# y = data["CC Share [%]"]
# ax.plot(
#     [x.min(), x.max()],
#     [y.min(), y.max()],
#     ls="--",
#     color="grey",
#     alpha=0.5,
#     zorder=-1,
# )
# ax.set_xlabel("Capacity Factor [%]")
# sns.despine()
# fig.savefig(snakemake.output[0], dpi=300)
