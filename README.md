
# maxRgain

<!-- badges: start -->
<!-- badges: end -->

maxRgain implements an integer programming-based method for optimizing
genetic gain in polyclonal selection, where the goal is to select a
group of genotypes that collectively meet multi-trait selection
criteria. The method relies on genetic effect predictors obtained from
mixed models and is demonstrated using grapevine breeding data, although
it is broadly applicable to other species and breeding contexts. This
package offers a ready-to-use implementation of the method proposed by
Surgy et al. (2025), allowing breeders to perform group-based
multi-trait selection in a consistent and computationally efficient
manner. It supports balanced genetic improvement across traits and
fosters transparency in breeding decision-making.

## Key Features

- Integer programming approach to group selection under multiple trait
  constraints

- Selection based on genetic effect predictors from mixed models

- Examples with real grapevine breeding data

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
#> $gain
#>         yd       pa       ta N.Clones
#> 1 24.95480 3.017470 3.049012       12
#> 2 25.49744 3.000773 3.011566       11
#> 3 26.57930 3.138235 3.068115       10
#> 4 27.74793 3.004891 3.103729        9
#> 5 27.85813 3.021805 3.101100        8
#> 6 28.80319 3.231580 3.120158        7
#> 
#> $selected
#>       12    11    10     9     8     7
#> 1  GV060 GV060 GV060 GV060 GV060 GV081
#> 2  GV079 GV080 GV081 GV079 GV080 GV088
#> 3  GV080 GV088 GV088 GV081 GV081 GV095
#> 4  GV081 GV089 GV093 GV088 GV088 GV097
#> 5  GV088 GV093 GV095 GV093 GV093 GV127
#> 6  GV089 GV097 GV097 GV097 GV097 GV140
#> 7  GV090 GV127 GV127 GV127 GV127 GV146
#> 8  GV097 GV133 GV133 GV140 GV146      
#> 9  GV127 GV140 GV140 GV146            
#> 10 GV133 GV144 GV146                  
#> 11 GV140 GV146                        
#> 12 GV146
```

For detailed guidance, see the package vignette:

``` r
vignette("maxRgain", package = "maxRgain")
```

## Citation

The method and the software can be cited as follows:

**Method:**  
Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful
tool for polyclonal selection in ancient grapevine varieties. Theor Appl
Genet 138, 122 (2025). <https://doi.org/10.1007/s00122-025-04885-0>

**Software (In preparation):**  
Surgy, S., Cadima, J. & Gonçalves, E. maxRgain - A package to perform
polyclonal selection. Manuscript in preparation.

## License

This package is released under the GPL-3 License.

## Contact and Bug reports

For suggestions, and bug reports use the link:

<https://github.com/soniasurgy/maxRgain/issues>
