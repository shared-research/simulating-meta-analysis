
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Understanding meta-analysis through data simulation with applications to power analysis

<!-- badges: start -->

[<img alt="alt_text" src="https://img.shields.io/badge/OSF-https%3A%2F%2Fosf.io%2F54djn%2F-blue" />](https://osf.io/54djn/)
<!-- badges: end -->

This repository contains the code to reproduce the analysis of the paper
**Understanding meta-analysis through data simulation with applications
to power analysis** by Gambarota Filippo and Gianmarco Altoè. The
project can also be found on the Open Science Framework
<https://osf.io/54djn/> with a preprint uploaded on PsyArXiv
<https://psyarxiv.com/br6vy/>.

## Usage

Functions used in the paper and supplementary materials are documented
and can be used opening the `simulating-meta-analysis.Rproj` file and
using the `devtools::load_all()` function (this requires the `devtools`
package installed). Then each function will be available as using a
standard R package.

## Documents

The repository (`documents/` folder) contains the paper and
supplementary materials created using R Markdown and the
[`papaja`](https://cran.r-project.org/web/packages/papaja/index.html)
package.

| **Document**            | **Rendered**                     | **Source**                       |
|-------------------------|----------------------------------|----------------------------------|
| Paper                   | [pdf](documents/paper/paper.pdf) | [rmd](documents/paper/paper.Rmd) |
| Supplementary Materials | [pdf](documents/paper/suppl.pdf) | [rmd](documents/paper/suppl.Rmd) |

## Repository

The project is organized as follows:

- `R`: contains the simulation and utility functions used in the paper
  and the supplementary materials
- `documents`: contains the `.Rmd` files to reproduce the paper and
  supplementary materials manuscripts
- `scripts`: contains the script used to produce complex objects and
  figures
  - `scripts\r-img.R`: create the paper figures. Figures that are not
    created in this script using `ggplot2` are created with the
    open-source program [*Inkscape*](https://inkscape.org/) and
    contained in the `../img/` folders.
- `objects`: contains the result of objects created by the `scripts/`
  scripts.

## Testing

The custom functions are tested within the `scripts/test.R` script. The
idea is to use the same simulation codes of the paper increasing the $k$
and $n$ values and check the recovery of simulated parameters using the
`check_sim()` function.

# Session

## Session info

    #> version: R version 4.3.1 (2023-06-16 ucrt)
    #> os: Windows 10 x64 (build 19045)
    #> system: x86_64, mingw32
    #> ui: RTerm
    #> language: (EN)
    #> collate: English_United States.utf8
    #> ctype: English_United States.utf8
    #> tz: Europe/Rome
    #> pandoc: 3.1.1 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

## Packages

- `dplyr` (1.1.2)
- `tidyr` (1.3.0)
- `ggplot2` (3.4.2)
- `metafor` (4.2.0)
- `MASS` (7.3.60)
- `cli` (3.6.1)
- `rmarkdown` (2.23)
- `papaja` (0.1.1)
- `knitr` (1.43)
- `devtools` (2.4.5)
- `here` (1.0.1)
- `kableExtra` (1.3.4.9000)
- `latex2exp` (0.9.6)
- `ltxplot` (1.0.0)
- `tibble` (3.2.1)
- `ggExtra` (0.10.0)
- `ggthemes` (4.2.4)
- `gridExtra` (2.3)
- `purrr` (1.0.2)
- `bookdown` (0.35)
- `broom` (1.0.5)
- `renv` (1.0.1)
- `stringr` (1.5.0)
- `simulatingMeta` (1.0.0)
- `sessioninfo` (1.2.2)
- `cowplot` (1.1.1)
- `distributional` (0.3.2)
- `ggdist` (3.3.0)
- `tidyverse` (2.0.0)
