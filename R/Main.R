#' Integer programming calculations
#'
#' @param traits A vector with the names of the columns in the data corresponding to the target traits to be optimized, i.e., those included in the objective function.
#' @param ref Name of the reference column (e.g., genotype ID). Defaults to the first column.
#' @param clmin An integer specifying the 	minimum group size.
#' @param clmax An integer specifying the maximum group size. If omitted, equal to `clmin`.
#' @param dmg A `data.frame` with three columns defining constraints: trait names; constraints signs (`">="`, `"<="` or `"=="`); and  _right-hand side_ values of the constraints.
#' @param meanvec Vector of overall means per trait; if omitted, data are assumed to be already normalized by the mean.
#' @param criteria Vector with the criterion for each trait: 1 to increase, -1 to decrease. If omitted, all traits are assumed to be increasing.
#' @param data  A data frame comprising the input data consisting of the Empirical Best Linear Unbiased Predictors (EBLUPs) of genotypic effects, which serve as the basis for the selection procedure.
#' @return
#' A list with the following components:
#'  - `gain`  with the gains of the several traits in each dimension
#'  - `selected`  with the reference os the clones selected in the group of each dimension
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). https://doi.org/10.1007/s00122-025-04885-0
#' @export
#' @examples
#' # The order of elements in the vectors `traits`, `meanvec`, and
#' #`criteria` does not affect the results.
#' # For example, `c("yd", "pa", "ta", "ph", "bw")` produces the same outcome as
#' # `c("pa", "yd", "ta", "ph", "bw")`only in a different order.
#' mymeanvec <- c(yd = 3.517, pa = 12.760, ta = 4.495, ph = 3.927, bw = 1.653)
#' mytraits <- c("yd", "pa",  "ta", "ph", "bw")
#' mydmg <- data.frame(
#'   lhs = c("yd", "pa", "ta", "ph", "bw"),
#'   rel = c(">=", ">=", ">=", ">=", ">="),
#'   rhs = c(20, 3, 3, 1, 2),
#'   stringsAsFactors = FALSE
#'   )
#' mycriteria <- c(yd = 1, pa = 1, ta = 1, ph = -1, bw = -1)
#' selections <- polyclonal(
#'    traits = mytraits,
#'    clmin = 7,
#'    clmax = 20,
#'    dmg = mydmg,
#'    meanvec = mymeanvec,
#'    criteria = mycriteria,
#'    data = Gouveio
#'    )
#' selections
polyclonal <- function(traits, ref = NULL, clmin = 7, clmax,  dmg = NULL, meanvec = NULL, criteria = NULL, data)
{
  if (length(traits) < 1) stop("There must be at least one trait in selection")
  if (!is.null(ref) && length(ref) != 1) stop("There must be only one reference column (ref).")
  if (is.null(ref)) ref <- names(data)[1]

  if (!all(traits %in% names(data))) stop("Some columns in 'traits' do not exist in data.")
  if (!all(ref %in% names(data))) stop("The reference column (ref) does not exist in data.")

  # Handle constraints
  const <- dmg[,1]
  relation <- dmg[,2]
  rhs <- dmg[,3]


  if (!all(const %in% names(data))) stop("Some constraint columns (dmg) do not exist in data.")

  # Trait + constraint columns
  cols <- unique(c(traits, const))


  # Check meanvec and criteria lengths
  if (!is.null(meanvec) && length(meanvec) != length(cols)) {
    stop("Length of 'meanvec' must match number of traits + constraints.")
  }
  if (!is.null(criteria) && length(criteria) != length(cols)) {
    stop("Length of 'criteria' must match number of traits + constraints.")
  }
  if (is.null(criteria)) criteria <- stats::setNames(rep(1, length(cols)), cols)
  if (is.null(meanvec)) meanvec <- stats::setNames(rep(1, length(cols)), cols)

  # Creating auxsum
  auxsum <- data.frame(col1=cols)

  if (!is.null(dmg)){
    idx <- match(auxsum$col1, dmg[[1]])
    auxsum$Crit <- criteria [idx]
    auxsum$Mean <- meanvec[idx]
    auxsum$DMG <- dmg[[3]][idx]
  }

  # Normalize data
  auxeblups <- data[, cols, drop = FALSE]
  auxeblups <- norm_eblup(auxeblups, cols, meanvec, criteria)
  auxeblups <- data.frame(data[,ref, drop = FALSE],auxeblups)
  if(missing(clmax)) clmax <- clmin
  if (clmax < clmin) stop("'clmax' must be greater than or equal to 'clmin'.")

  clmin <- as.integer(clmin)
  clmax <- as.integer(clmax)
  if (clmin < 1) stop("'clmax' must be greater than or equal to 'clmin'.")

  return(ipp(traits=traits, ref = ref, clmin, clmax, const = const,  relation = relation, rhs = rhs, dmg = dmg, auxsum = auxsum, data = auxeblups ))
}
