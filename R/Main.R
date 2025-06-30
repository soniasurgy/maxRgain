#' Integer programming calculations
#'
#' @param traits A vector with the names of the columns in the data corresponding to the target traits to be optimized, i.e., those included in the objective function.
#' @param ref A vector with the name of the column containing the reference identifiers of the clones.
#' @param clmin An integer specifying the number of clones to include in the smallest polyclonal group.
#' @param clmax An integer specifying the number of clones to include in the largest polyclonal group.
#' @param dmg Desired minimum gain for the traits in the constraints.
#' @param meanvec A vector with the phenotypic mean values, with names corresponding to column names (e.g., trait1 = mean1, trait2 = mean2). If not provided, it is assumed the values are already divided by the mean.
#' @param criteria A vector with the selection criteria, with values of 1 or -1 depending on whether an increase or decrease in the trait is desired, respectively (e.g., trait1 = 1, trait2 = -1). If not provided, it is assumed that an increase is wanted.
#' @param data  A data frame comprising the input data consisting of the Estimated Best Linear Unbiased Predictors (EBLUPs) of genotypic effects, which serve as the basis for the selection procedure.
#' @return list with 2 objects
#'     $gain  With the gains of the several traits in each dimension
#'     $selected  with the reference os the clones selected in the group of each dimension
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
#' selections$gain
#' selections$selected
polyclonal <- function(traits, ref = NULL, clmin = 7, clmax,  dmg = NULL, meanvec = NULL, criteria = NULL, data)
{
  if (length(traits) < 1){
    stop("There must be at least one trait in selection")
  }
  if (!is.null(ref)){
    if (length(ref) != 1){
      stop("There must be only one ref")
    }
  }
  const <- dmg[,1]
  relation <- dmg[,2]
  rhs <- dmg[,3]
  if (is.null(ref)){
    ref <- names(data)[1]
  }
  if (!all(traits %in% names(data))) stop("Some columns in traits don't exist.")
  if (!all(const %in% names(data))) stop("Some columns in const don't exist.")
  if (!all(ref %in% names(data))) stop("The column in ref doesn't exist.")
  cols <- unique(c(traits, const))
  if (!is.null(meanvec)){
    if (length(meanvec) != length(cols)){
      stop("Wrong length of vector meanvec")
    }
  }
  if (!is.null(criteria)){
    if (length(criteria) != length(cols)){
      stop("Wrong length of vector criteria")
    }
  }
  auxeblups <- data[, cols, drop = FALSE]
  auxeblups <- norm_eblup(auxeblups, cols, meanvec, criteria)
  auxeblups <- data.frame(data[,ref, drop = FALSE],auxeblups)
  if(missing(clmax)){
    clmax <- clmin
  }
  if (clmax < clmin){
    stop("clmax must be grater then clmin")
  }
  intclmin <- utils::type.convert(clmin, as.is=integer)
  intclmax <- utils::type.convert(clmax, as.is=integer)
  if (intclmin < 1){
    stop("clmin must be a positive integer")
  }
  return(ipp(traits=traits, ref = ref, clmin, clmax, const = const,  relation = relation, rhs = rhs, data = auxeblups))
}
