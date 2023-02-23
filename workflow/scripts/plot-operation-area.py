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
        "plot_operation_area",
        simpl="",
        lv=1.2,
        clusters=181,
        opts="",
        sector_opts="Co2L0-365H-T-H-B-I-A-solar+p3-linemaxext15-seq200",
        planning_horizons=2050,
        kind="electricity",
    )


kind = snakemake.wildcards.kind
config = snakemake.config
which = "operation"


df = {}
for path in snakemake.input.networks:
    n = import_network(path)
    prod = get_carrier_production(n, kind, config, which).groupby(level=1).sum()
    cons = get_carrier_consumption(n, kind, config, which).groupby(level=1).sum()

    if (prod.sum() / cons.sum()).round(3) != 1:
        print(f"Warning: {kind} production and consumption are not equal.")

    sec_opts = path.split("/")[-1].split("_")[-2].split("-")
    sequestration = int([s for s in sec_opts if "seq" in s][0][3:])
    if "CF+sector+co2network+false" in sec_opts:
        co2_network = "Without Carbon Network"
    else:
        co2_network = "With Carbon Network"

    key = (co2_network, sequestration)

    df[key] = pd.concat([prod, cons], keys=["Production", "Consumption"])
df = pd.concat(df, axis=1)


norm = 1e6
unit = "TWh" if kind not in ["carbon", "co2"] else "Mt"
nice_name = n.carriers.nice_name
colors = n.carriers.color.dropna().rename(nice_name)
sort_by_color = (
    lambda df: df.assign(color=colors[df.index])
    .sort_values(by="color")
    .drop("color", axis=1)
)

keys = df.columns.unique(level=0)
for k, output in zip(keys, snakemake.output):

    production = (
        df.loc["Production", k].rename(index=nice_name).div(norm).pipe(sort_by_color).T
    )
    consumption = (
        df.loc["Consumption", k].rename(index=nice_name).div(norm).pipe(sort_by_color).T
    )

    production = production.loc[:, production.sum() > 0.1]
    consumption = consumption.loc[:, consumption.sum() > 0.1]

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    production.plot(kind="area", ax=ax, color=colors.to_dict(), alpha=0.8)
    consumption.mul(-1).plot(kind="area", ax=ax, color=colors.to_dict(), alpha=0.8)
    ax.axhline(0, color="k", lw=1)
    ax.set_xlim(production.index.min(), production.index.max())
    ax.set_ylim(-consumption.sum(1).max() * 1.1, production.sum(1).max() * 1.1)
    ax.set_xlabel("Sequestration Potential [Mt]")
    ax.set_title(f"{kind.title()} Balance {k}")
    ax.grid(axis="y", alpha=0.5)

    # split labels and handles into production and consumption
    handles, labels = ax.get_legend_handles_labels()
    pcarriers = production.columns
    ccarriers = consumption.columns
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
        bbox_to_anchor=(1, 1),
        frameon=False,
        ncol=2,
        title="Production",
        labelcolor="k",
    )
    fig.add_artist(legend)

    legend = ax.legend(
        chandles[::-1],
        clabels[::-1],
        loc="lower left",
        bbox_to_anchor=(1, 0),
        frameon=False,
        ncol=2,
        title="Consumption",
        labelcolor="k",
    )
    fig.add_artist(legend)

    sns.despine()
    ax.set_ylabel(unit)

    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
