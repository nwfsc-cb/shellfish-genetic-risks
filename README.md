
<!-- README.md is generated from README.Rmd. Please edit that file -->

# shellfishrisks

<!-- badges: start -->
<!-- badges: end -->

Repo for hosting R package for shellfish genetic risk modeling

This package allows for agent based modeling of shellfish genetic risk
modeling. The model itself runs in Python, but users can use this
package to both run and process the results in R.

## Installation

To install shellfishrisks from [GitHub](https://github.com/) run

``` r
# install.packages("devtools")
devtools::install_github("nwfsc-cb/shellfish-genetic-risks")
```

Note that this requires installing the `devtools` package prior to
installing `shellfishrisks`

This package also requires a working Python 3.X installation. The model
is currently built assuming users have an Anaconda installation of
Python. Users unfamiliar with installing Python should install Anaconda
and Python 3.X following the instructions
[here](https://docs.anaconda.com/anaconda/install/). If you are familiar
enough with Python to have strong opinions about the type of Python
installation you would like (i.e. something other than Anaconda) we
shall assume you can do so without our guidance.

We won’t kid you. Users unfamiliar with Python installation may find the
process extremely daunting . We have tried to make this process as
painless as possible. The package is designed assuming use of Anaconda.
We then interface with Python through R using the `reticulate` package.
By default, `reticulate` automatically creates and uses a Conda
(Anaconda) environment called `r-reticulate`. The packages attempts to
install any missing and required Python packages the first time the
package is run, so if you have not previously installed the required
Python packages you will need a working internet connection.

See the example below for the workflow needed to use this package

## Example

This example shows you the sequence of steps needed to run
`shellfishrisks`.

``` r
# For various reasons, if think you need to install Python libraries, you need to run this BEFORE calling library(shellfishrisks)

shellfishrisks::install_pylibs()

library(shellfishrisks) # now load shellfishrisks after installing required python libraries

# Set options
batch_name <- "dev"

reps <- 1

coreid <- 1

# All the _years should ideally be set to 50. However, that takes >24 hours to run. For testing, we recommend setting years to 1 to ensure the installation is working properly. 
pre_farm_years <- 1 

farm_years <- 1

post_farm_years <- 1

# run shellfishrisk model, should take about 10 minutes with these settings
eggs(
  batch = batch_name,
  reps = as.integer(reps),
  coreid = as.integer(coreid),
  pre_farm_years = as.integer(pre_farm_years),
  farm_years = as.integer(farm_years),
  post_farm_years = as.integer(post_farm_years)
)

clean_shellfish(batch = batch_name) # move the results to a folder of the form results/{batch}

results <- process_shellfish(batch = batch_name, results_dir = file.path("results",batch_name)) # read the results stored in .txt files into a list object

str(results)
```
