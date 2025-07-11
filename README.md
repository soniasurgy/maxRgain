
# maxRgain

maxRgain implements an Integer Programming-based method for optimizing
genetic gain in polyclonal selection, where the goal is to select a
group of genotypes that jointly meet multi-trait selection criteria. The
method uses predictors of genotypic effects obtained from the fitting of
mixed models. Its application is demonstrated with grapevine data, but
is applicable to other species and breeding contexts.

## Key Features

- Integer programming approach to group selection under multiple trait
  constraints

- Selection based on predictors of genotypic effects from mixed models

- Examples with real data from field trials of grapevine

- The method can be generalized to a variety of crops and breeding
  programs

- Includes example dataset and reproducible workflows

## Installation

To install the development version of maxRgain from
[GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("soniasurgy/maxRgain")
```

To install from CRAN:

``` r
install.packages("maxRgain")
```

## Example

This is a basic example which shows how to solve a common problem:

``` r
library(maxRgain)
data("Gouveio")
mydmg <- data.frame(
  lhs = c("yd", "pa", "ta"),
  rel = c(">=", ">=", ">="),
  rhs = c(20, 3, 3),
  stringsAsFactors = FALSE
  )

polyclonal(
  traits = c("yd", "pa", "ta"),
  clmin = 7,
  clmax = 12,
  dmg = mydmg,
  meanvec = c(yd = 3.517, pa = 12.760, ta = 4.495),
  data = Gouveio
  )
#> Predited genetic gains as a % of the overall mean
#> 
#> $gain 
#>  Group.Size        yd       pa       ta
#>          12  87.76604 38.50292 13.70531
#>          11  89.67451 38.28986 13.53699
#>          10  93.47938 40.04388 13.79118
#>           9  97.58947 38.34240 13.95126
#>           8  97.97704 38.55823 13.93945
#>           7 101.30081 41.23496 14.02511
#> 
#> 
#> Selected genotypes (per group size)
#> 
#> $selected
#>     12    11    10     9     8     7
#>  GV060 GV060 GV060 GV060 GV060 GV081
#>  GV079 GV080 GV081 GV079 GV080 GV088
#>  GV080 GV088 GV088 GV081 GV081 GV095
#>  GV081 GV089 GV093 GV088 GV088 GV097
#>  GV088 GV093 GV095 GV093 GV093 GV127
#>  GV089 GV097 GV097 GV097 GV097 GV140
#>  GV090 GV127 GV127 GV127 GV127 GV146
#>  GV097 GV133 GV133 GV140 GV146      
#>  GV127 GV140 GV140 GV146            
#>  GV133 GV144 GV146                  
#>  GV140 GV146                        
#>  GV146
```

For detailed guidance, see the package vignette:

``` r
vignette("maxRgain", package = "maxRgain")
```

## Citation

The underlying method and this package can be cited as follows:

**Method:**  
Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful
tool for polyclonal selection in ancient grapevine varieties. Theor Appl
Genet 138, 122 (2025). <https://doi.org/10.1007/s00122-025-04885-0>

**Package (In preparation):**  
Surgy, S., Cadima, J. & Gonçalves, E. maxRgain - A package to maximize
genetic gains of polyclonal selection. Manuscript in preparation.

## License

This package is released under the GPL-3 License.

## Contact and Bug reports

For suggestions, and bug reports use the link:

<https://github.com/soniasurgy/maxRgain/issues>
