
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Understanding meta-analysis through data simulation, with applications to power analysis

<!-- badges: start -->

[<img alt="alt_text" src="https://img.shields.io/badge/OSF-https%3A%2F%2Fosf.io%2F54djn%2F-blue" />](https://osf.io/54djn/)
<!-- badges: end -->

This repository contains the code to reproduce the analysis of the paper
**Understanding meta-analysis through data simulation, with applications
to power analysis** by Gambarota Filippo and Gianmarco Altoè. The
project can also be found on the Open Science Framework
<https://osf.io/54djn/>.

## Repository

The project is organized as follows:

- `R`: contains the simulation functions used in the paper and
  supplementary materials
- `docs`: contains the `.Rmd` files to reproduce the paper and
  supplementary materials manuscripts
- `scripts`: contains the script used to produce complex objects and
  figures
  - `scripts\r-img.R`: create the paper figures. Figures that are not
    produced in this script using `ggplot2` are created with the
    open-source program [*Inkscape*](https://inkscape.org/) and
    contained in the `../img/` folders.
- `objects`: contains the result of objects created by the `scripts/`
  scripts.

The project is organized and structured as an R package, thus after
opening the `simulating-meta-analysis.Rproj` file and installing the
`devtools` package, calling `devtools::load_all()` will load each
function within the `R/` folder.
