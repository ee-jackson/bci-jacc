# bci-jacc

This repository contains the [research compendium](https://research-compendium.science) for our in-prep paper: 

Eleanor E. Jackson, Patrick G. Cannon, Matthew P. Greenwell, Sofia Gripenberg. **Density-dependent pre-dispersal seed predation and fruit set in a tropical tree: spatial scaling of plant-insect interactions.** 

## Contents:

### [`code/`](code/)
The [`code/`](code/) directory contains two subdirectories:

[`exploration/`](code/exploration/) contains R Markdown notebook files that are used for exploratory analysis. They produce figures and explanatory documents.  All `.Rmd` files in [`exploration/`](code/exploration/) are knitted to `.md` `github_documents` to make the GitHub repo browsable.

[`scripts/`](code/scripts/) contains action scripts, i.e. all the code for cleaning, combining, and analysing the data. Each `.R` script has a summary at the top of what it does. The scripts are numbered in the order in which they would typically be run.

### `data/`
The original data is stored in the `data/raw/` subdirectory. Any data that is produced using code is stored in `data/clean/`. Data will be made available on publication of the research paper.

### [`output/`](output/)
The [`output/`](output/) directory contains the subdirectories [`figures/`](output/figures/) and [`results/`](output/results/), which contain the figures used in the paper and other output from analyses, respectively.

### [`docs/`](docs/)
The [`docs/`](docs/) directory contains the [data dictionary](docs/data-dictionary.md) (i.e. metadata), [project notebook](docs/project-notebook.md), [protocols](docs/protocols.md) and any other relevant documents.
