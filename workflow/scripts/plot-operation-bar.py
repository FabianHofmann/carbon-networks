import os
from common import (
    import_network,
    mock_snakemake,
    get_carrier_consumption,
    get_carrier_production,
)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="white", context="paper")


if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_operation_bar",
        simpl="",
        lv=1.2,
        clusters=181,
        opts="",
        sector_opts="Co2L0-365H-T-H-B-I-A-solar+p3-linemaxext15-seq200",
        planning_horizons=2050,
    )


n = import_network(snakemake.input.network)
kinds = snakemake.config["constants"]["kinds"]
config = snakemake.config
which = "operation"


df = {}
for kind in kinds:
    prod = get_carrier_production(n, kind, config, which).groupby(level=1).sum()
    cons = get_carrier_consumption(n, kind, config, which).groupby(level=1).sum()

    if (prod.sum() / cons.sum()).round(3) != 1:
        print(f"Warning: {kind} production and consumption are not equal.")

    df[kind] = pd.concat([prod, cons], keys=["Production", "Consumption"])
df = pd.concat(df)


fig, axes = plt.subplots(len(kinds), 1, figsize=(14, 14), sharex=True)
colors = n.carriers.color.dropna()
nice_name = n.carriers.nice_name

for kind, ax in zip(kinds, axes):

    norm = 1e6
    unit = "TWh" if kind not in ["carbon", "co2"] else "Mt"

    sdf = df.div(norm).loc[kind].unstack(level=1).loc[["Production", "Consumption"]]
    sdf.loc["color"] = colors[sdf.columns]
    sdf = sdf.sort_values(axis=1, by="color").drop("color", axis=0)
    sdf = sdf.rename(columns=nice_name)
    sdf.plot.bar(
        stacked=True,
        ax=ax,
        legend=False,
        rot=0,
        color=colors.rename(nice_name).to_dict(),
        alpha=0.8,
    )
    ax.set_title(kind.title())
    ax.grid(axis="y", alpha=0.5)

    # split labels and handles into production and consumption
    handles, labels = ax.get_legend_handles_labels()
    pcarriers = sdf.loc["Production"].dropna()[lambda ds: ds >= 0.1].index
    ccarriers = sdf.loc["Consumption"].dropna()[lambda ds: ds >= 0.1].index
    plabels, phandles = [], []
    clabels, chandles = [], []
    for label, handle in zip(labels[::-1], handles[::-1]):
        if label in pcarriers:
            plabels.append(label)
            phandles.append(handle)
        if label in ccarriers:
            clabels.append(label)
            chandles.append(handle)

    legend = ax.legend(
        phandles,
        plabels,
        loc="upper left",
        bbox_to_anchor=(1, 1.1),
        frameon=False,
        ncol=1,
        title="Production",
        labelcolor="k",
    )
    fig.add_artist(legend)

    legend = ax.legend(
        chandles,
        clabels,
        loc="upper left",
        bbox_to_anchor=(1.7, 1.1),
        frameon=False,
        ncol=1,
        title="Consumption",
        labelcolor="k",
    )
    fig.add_artist(legend)

    sns.despine()
    ax.set_ylabel(unit)

fig.tight_layout()
fig.savefig(snakemake.output[0], dpi=300, bbox_inches="tight")
