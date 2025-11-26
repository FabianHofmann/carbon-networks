# Combining the use of CO₂ and H₂ networks benefits carbon management in Europe

This repository contains the official workflow to reproduce the results of the paper:

> **Combining the use of CO₂ and H₂ networks benefits carbon management in Europe**
> Fabian Hofmann, Christoph Tries, Fabian Neumann, Elisabeth Zeyen, Tom Brown
> *Nature Energy* (2025)
> DOI: [10.1038/s41560-025-01753-5](https://doi.org/10.1038/s41560-025-01753-5)
> Preprint: [arXiv:2402.19042](https://arxiv.org/abs/2402.19042)

The study investigates how hydrogen and carbon dioxide transport networks can support a climate-neutral European energy system. It finds that a hydrogen network is more cost-effective for transporting hydrogen to demand centers, while both networks complement each other when deployed together—the CO₂ network encourages distributed biomass capture and reduces dependence on direct air capture.

## Funding

This project was funded by [Breakthrough Energy](https://www.breakthroughenergy.org/).

## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

## Run the analysis

    snakemake -call

This will run all analysis steps to reproduce results and eventually build the report.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake -c1 --use-conda -f dag

## Repo structure

* `config`: configurations used in the study
* `data`: place for raw data
* `report`: contains all files necessary to build the report; plots and result files are generated automatically
* `workflow`: contains the Snakemake workflow
* `build`: will contain all results (does not exist initially)

## Citation

If you use this code, please cite:

```bibtex
@article{hofmann2025co2h2,
  title={Combining the use of CO2 and H2 networks benefits carbon management in Europe},
  author={Hofmann, Fabian and Tries, Christoph and Neumann, Fabian and Zeyen, Elisabeth and Brown, Tom},
  journal={Nature Energy},
  year={2025},
  doi={10.1038/s41560-025-01753-5}
}
```

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
