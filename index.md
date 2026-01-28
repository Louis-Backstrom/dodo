# Dodo

This package includes many extinction models published in the ecology
and palaeontology literature. It is still in active development.

## Installing

    ## install and load devtools (if needed)
    install.packages("devtools") # if not already installed
    library(devtools)

    ## install dodo from GitHub and load
    install_github("Louis-Backstrom/dodo")
    library(dodo)

## Getting Started

    ## access list of package functions
    ?dodo

    ## explore monk seal sightings dataset from Solow 1993a
    monk_seal

    ## fit the monk seal model from Solow 1993a
    SO93F1(monk_seal, alpha = 0.05, init.time = 1915, test.time = 1992)

## Dependencies

Certain functions in this package require
[JAGS](https://mcmc-jags.sourceforge.io/) (accessed via `rjags`). This
should be installed prior to running any of the functions which
incorporate Gibbs sampling.

## Contact

Louis Backstrom (<ljb38@st-andrews.ac.uk>)

Package website (work in progress):
[louis-backstrom.github.io/dodo/](https://louis-backstrom.github.io/dodo/louis-backstrom.github.io/dodo/)
