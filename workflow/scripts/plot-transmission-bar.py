import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from common import (
    get_transmission_links,
    import_network,
    mock_snakemake,
    sort_rows_by_diff,
)

alpha = 1
region_alpha = 0.8

if os.path.dirname(os.path.abspath(__file__)) == os.getcwd():
    snakemake = mock_snakemake(
        "plot_transmission_bar", ext="pdf", clusters=90, comparison="default"
    )

sns.set_theme(**snakemake.params["theme"])

df = {}
objectives = {}
for path in snakemake.input.networks:
    n = import_network(path)

    is_transport = get_transmission_links(n)
    transport_carriers = [
        *n.links.carrier[is_transport].unique(),
        *n.lines.carrier.unique(),
    ]
    transport_carriers = n.carriers.nice_name[transport_carriers]

    transmission = n.statistics.supply(comps=n.branch_components)
    transmission = transmission.reindex(transport_carriers, level=1)
    transmission = transmission.droplevel(0)

    key = snakemake.params.labels[n.meta["wildcards"]["run"]]

    df[key] = transmission

data = pd.concat(df, axis=1)
colors = n.carriers.set_index("nice_name").color.to_dict()

norm = 1e9
data = sort_rows_by_diff(data).div(norm)

fig, ax = plt.subplots(
    1, 1, figsize=snakemake.params.settings["figsize"], layout="constrained"
)


defaults = dict(kind="bar", stacked=True, lw=0, rot=0, alpha=0.8)

kwargs = {**defaults, **snakemake.params.settings.get("kwargs", {})}
data[::-1].T.plot(ax=ax, color=colors, **kwargs)

ax.axhline(0, color="k", lw=1)
ax.set_ylabel("Transport Volume [pWh] / [Mt]")
ax.grid(axis="y", alpha=0.5)
handles, labels = ax.get_legend_handles_labels()
ax.legend().remove()


sns.despine()
fig.legend(
    handles[::-1],
    labels[::-1],
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    frameon=False,
    ncol=1,
)
fig.savefig(snakemake.output[0], bbox_inches="tight", dpi=300)
