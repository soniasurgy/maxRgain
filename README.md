
# maxRgain

maxRgain implements an Integer Programming-based method for optimizing
genetic gain in polyclonal selection, where the goal is to select a
group of genotypes that jointly meet multi-trait selection criteria.
Polyclonal selection, as opposed to the selection of a single genotype,
helps to preserve genetic diversity and provides greater environmental
stability, while still obtaining high genetic gains. The method uses
predictors of genotypic effects obtained from the fitting of mixed
models. Its application is demonstrated with grapevine data, but is
applicable to other species and breeding contexts.

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

## Dependencies

This package depends on the following R packages:

- **lpSolve**: for solving integer linear programming problems.
- **stats** and **utils**: base R packages used for common statistical
  and utility functions.

## Example

In this example, we wish to improve yields by at least 20%, and the
potential alcohol and total acidity by at least 3%. We wish to obtain
group sizes between 7 and 12 genotypes.

``` r
library(maxRgain)

polyclonal(
  traits = c("yd", "pa", "ta"),
  clmin = 7,
  clmax = 12,
  dmg = data.frame(
    lhs = c("yd", "pa", "ta"),
    rel = c(">=", ">=", ">="),
    rhs = c(20, 3, 3)
    ),
  meanvec = c(yd = 3.517, pa = 12.760, ta = 4.495),
  data = Gouveio
  )
#> Predited genetic gains as a % of the overall mean
#> 
#> $gain 
#>  Group.Size       yd       pa       ta
#>          12 24.95480 3.017470 3.049012
#>          11 25.49744 3.000773 3.011566
#>          10 26.57930 3.138235 3.068115
#>           9 27.74793 3.004891 3.103729
#>           8 27.85813 3.021805 3.101100
#>           7 28.80319 3.231580 3.120158
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
