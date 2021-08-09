
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

This example shows you the sequence of steps needed to install requirede
Python packages and run `shellfishrisks`.

``` r
set.seed(42)
# For various reasons, if think you need to install Python libraries, you need to run this BEFORE calling library(shellfishrisks)

shellfishrisks::install_pypkgs() # installs any missing and required Python packages
#> shellfishrisks is using Conda environment 'r-reticulate'

library(shellfishrisks) # now load shellfishrisks after installing required python packages
# If you get an error about a missing Python package, restart R, and then run shellfishrisks::install_pypgks() BEFORE running library(shellfishrisks)

library(ggplot2)

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

# Set options
batch_name <- "dev"

reps <- 1

coreid <- 1

# All the _years should ideally be set to 50. However, that takes >24 hours to run. For testing, we recommend setting years to 1 to ensure the installation is working properly. 
pre_farm_years <- 1 

farm_years <- 1

post_farm_years <- 1


eggs(
  batch = batch_name,
  reps = reps,
  coreid = coreid,
  pre_farm_years = pre_farm_years,
  farm_years = farm_years,
  post_farm_years = post_farm_years,
  wild_N_init = 100
)

clean_shellfish(batch = batch_name) # move the results to a folder of the form results/{batch}
#> results moved to results/dev

results <- process_shellfish(batch = batch_name, results_dir = file.path("results",batch_name)) # read the results stored in .txt files into a list object

str(results)
#> List of 8
#>  $ AFs             :'data.frame':    500 obs. of  10 variables:
#>   ..$ coreid         : int [1:500] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Srep           : int [1:500] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Rep            : int [1:500] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Year           : int [1:500] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Subpop         : chr [1:500] "wild1" "wild1" "wild1" "wild1" ...
#>   ..$ Locus_index    : int [1:500] 0 0 1 1 2 2 3 3 4 4 ...
#>   ..$ Allele         : int [1:500] 0 1 0 1 0 1 0 1 0 1 ...
#>   ..$ Adaptive       : chr [1:500] "True" "True" "True" "True" ...
#>   ..$ Adv            : chr [1:500] "True" "False" "True" "False" ...
#>   ..$ AlleleFrequency: num [1:500] 0.9048 0.0952 0.8563 0.1437 0.9545 ...
#>  $ bstock_report   :'data.frame':    8 obs. of  4 variables:
#>   ..$ coreid            : int [1:8] 1 1 1 1 1 1 1 1
#>   ..$ Rep               : int [1:8] 0 0 0 0 0 0 0 0
#>   ..$ Year              : int [1:8] 1 1 1 1 1 1 1 1
#>   ..$ IndID_survivedAgeX: num [1:8] 68 96 42 89 65 2 48 64
#>  $ harvest         :'data.frame':    3 obs. of  5 variables:
#>   ..$ coreid       : int [1:3] 1 1 1
#>   ..$ Srep         : int [1:3] 1 1 1
#>   ..$ Rep          : int [1:3] 0 0 0
#>   ..$ Year         : int [1:3] 0 1 2
#>   ..$ IndsHarvested: int [1:3] 0 0 0
#>  $ life_hist_report:'data.frame':    124701 obs. of  12 variables:
#>   ..$ coreid       : int [1:124701] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Rep          : int [1:124701] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Year         : int [1:124701] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Month        : chr [1:124701] "Jan" "Jan" "Jan" "Jan" ...
#>   ..$ Subpop       : chr [1:124701] "wild1" "wild1" "wild1" "wild1" ...
#>   ..$ IndID        : num [1:124701] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..$ Sex          : int [1:124701] 1 2 2 2 1 2 2 2 1 2 ...
#>   ..$ CohortYear   : num [1:124701] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Age          : num [1:124701] 7.25 4.33 1.17 4.42 4.58 ...
#>   ..$ Mother_id    : num [1:124701] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Father_id    : num [1:124701] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ meanParentAge: num [1:124701] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ log             :'data.frame':    3 obs. of  9 variables:
#>   ..$ coreid    : int [1:3] 1 1 1
#>   ..$ rep       : int [1:3] 0 0 0
#>   ..$ farm_phase: chr [1:3] "pre-farm" "during-farm" "post-farm"
#>   ..$ year      : int [1:3] 0 1 2
#>   ..$ wild1_size: int [1:3] 725 1212 1614
#>   ..$ wild2_size: int [1:3] 701 1223 1648
#>   ..$ wild3_size: int [1:3] 726 1224 1619
#>   ..$ farm_size : int [1:3] NA 497 NA
#>   ..$ time_stamp: chr [1:3] "2021-08-06_09:37:41" "2021-08-06_09:40:14" "2021-08-06_09:43:31"
#>  $ pop_pair_rvars  :'data.frame':    36 obs. of  7 variables:
#>   ..$ coreid  : int [1:36] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Srep    : int [1:36] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Rep     : int [1:36] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Year    : int [1:36] 0 0 0 0 0 0 0 0 0 1 ...
#>   ..$ Pop_pair: chr [1:36] "wild1_wild2" "wild1_wild2" "wild1_wild2" "wild2_wild3" ...
#>   ..$ Rvar    : chr [1:36] "Fst" "Fst_aL" "Fst_nL" "Fst" ...
#>   ..$ Value   : num [1:36] 0.00947 0.00774 0.01039 0.00604 0.00434 ...
#>  $ pop_rvars       :'data.frame':    40 obs. of  7 variables:
#>   ..$ coreid: int [1:40] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Srep  : int [1:40] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ Rep   : int [1:40] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Year  : int [1:40] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Subpop: chr [1:40] "wild1" "wild1" "wild1" "wild1" ...
#>   ..$ Rvar  : chr [1:40] "mfit" "popsize" "het" "ar" ...
#>   ..$ Value : num [1:40] 0.996 835 1 1740 0.997 ...
#>  $ temp_rvars      :'data.frame':    3 obs. of  6 variables:
#>   ..$ coreid: int [1:3] 1 1 1
#>   ..$ Srep  : int [1:3] 1 1 1
#>   ..$ Rep   : int [1:3] 0 0 0
#>   ..$ Subpop: chr [1:3] "wild1" "wild2" "wild3"
#>   ..$ Rvar  : chr [1:3] "Fst" "Fst" "Fst"
#>   ..$ Value : num [1:3] 0.000528 0.000238 0.00018

results$AFs %>% 
  ggplot(aes(Year, AlleleFrequency, color = Locus_index, group = Locus_index)) + 
  geom_line() + 
  facet_grid(Subpop ~ Allele, labeller = label_both)
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
#> geom_path: Each group consists of only one observation. Do you need to adjust
#> the group aesthetic?
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r


results$pop_pair_rvars %>% 
  ggplot(aes(Year, Value, color = Pop_pair)) + 
  geom_line() + 
  facet_wrap(~ Rvar, labeller = label_both, scales = "free_y")
```

<img src="man/figures/README-example-2.png" width="100%" />
