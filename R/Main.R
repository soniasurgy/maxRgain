#' Integer programming calculations
#'
#' This function maximizes the predicted genetic gain in the selection of groups of genotypes based on the predictores of genotypic effects.
#'
#' @param traits A vector with the names of the columns in the data corresponding to the target traits to be optimized, i.e., those included in the objective function.
#' @param ref Name of the reference column (e.g., genotype ID). Defaults to the first column.
#' @param clmin An integer specifying the 	minimum group size.
#' @param clmax An integer specifying the maximum group size. If omitted, equal to `clmin`.
#' @param dmg A `data.frame` with three columns defining constraints: trait names; constraints signs (`">="`, `"<="` or `"=="`); and  _right-hand side_ values of the constraints.
#' @param meanvec A named numeric vector of overall means per trait; if omitted, data are assumed to be already normalized by the mean.
#' @param criteria A named numeric vector indicating the selection criterion for each trait: 1 for traits to be increased, -1 for traits to be decreased. If omitted, all traits are assumed to be selected for increase.
#' @param data  A data frame comprising the input data consisting of the Empirical Best Linear Unbiased Predictors (EBLUPs) of genotypic effects, which serve as the basis for the selection procedure.
#' @note The order of traits must be consistent across `traits`, `dmg`, `meanvec`, and `criteria`. Both `meanvec` and `criteria` must include values for all traits specified in `traits` and `dmg`.
#' @return
#' A list with the following components:
#'  - `gain`  with the gains of the several traits in each dimension
#'  - `selected`  with the reference os the clones selected in the group of each dimension
#'
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). https://doi.org/10.1007/s00122-025-04885-0
#' @export
#' @examples
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
  # Input checks
  if (length(traits) < 1) stop("There must be at least one trait in selection")
  if (!is.null(ref) && length(ref) != 1) stop("Only one reference column is allowed (ref).")
  if (is.null(ref)) ref <- names(data)[1]

  if (!all(traits %in% names(data))) stop("Some columns in 'traits' do not exist in `data`.")
  if (!all(ref %in% names(data))) stop("The reference column (ref) does not exist in `data`.")

  # Constraint setup
  const <- if (!is.null(dmg)) dmg[[1]] else NULL
  relation <- if (!is.null(dmg)) dmg[[2]] else NULL
  rhs <- if (!is.null(dmg)) dmg[[3]] else NULL


  if (!is.null(const) && !all(const %in% names(data))) {
    stop("Some constraint in `dmg` do not exist in 'data'.")
  }

  # Trait + constraint columns
  cols <- unique(c(traits, const))


  # Validate vector lengths
  if (!is.null(meanvec) && length(meanvec) != length(cols)) {
    stop("Length of 'meanvec' must contain values for all traits in `traits` and `constraints`.")
  }
  if (!is.null(criteria) && length(criteria) != length(cols)) {
    stop("Length of 'meanvec' must contain values for all traits in `traits` and `constraints`.")
  }

  # Set defaults
  if (is.null(criteria)) criteria <- stats::setNames(rep(1, length(cols)), cols)
  if (is.null(meanvec)) meanvec <- stats::setNames(rep(1, length(cols)), cols)

  # Create summary matrix
  auxsum <- data.frame(col1=cols)

  if (!is.null(dmg)){
    idx <- match(auxsum$col1, dmg[[1]])
    auxsum$Crit <- criteria [idx]
    auxsum$Mean <- meanvec[idx]
    auxsum$DMG <- rhs[idx]
  }

  # Normalize trait and constraint values
  auxeblups <- norm_eblup(data[, cols, drop = FALSE], cols, meanvec, criteria)
  auxeblups <- data.frame(data[, ref, drop = FALSE], auxeblups)

  # Check clone group size limits
  if(missing(clmax)) clmax <- clmin
  if (clmax < clmin || clmin < 1) stop("'clmax' must be >= 'clmin' and both >= 1.")

  clmin <- as.integer(clmin)
  clmax <- as.integer(clmax)

  # Objective function (sum of traits to maximize)
  objdir <- "max"
  auxobj <- auxeblups[, traits, drop = FALSE]
  fobj=0
  fobj <- Reduce(`+`, auxobj)

  # Constraint matrix
  if (is.null(const)) {
    constvalue <- data.frame(rep(1, nrow(auxeblups)))
    colnames(constvalue) <- "var"
    dirvalue <- "=="
  } else {
    constvalue <- data.frame(rep(1, nrow(auxeblups)), auxeblups[, const, drop = FALSE])
    colnames(constvalue) <- c("var", const)
    dirvalue <- c("==", relation)
  }

  # Vars for results
  clsel <- 0
  clselmax <- 0
  indsol <- 0
  gainout <- NULL
  auxeblups <- data.frame(data)

  # Iterate over group sizes
  for (i in clmax:clmin){
    rhsvalue  <- i
    if (!is.null(rhs)){
      for (a in 1:length(rhs)){
        rhsvalue <- c(rhsvalue,  rhs[a]*i/100)
      }
    }

    # Solver
    prob <- lpSolve::lp(
      direction = objdir,
      objective.in = fobj,
      const.mat = constvalue,
      const.dir=dirvalue,
      const.rhs=rhsvalue,
      transpose.constraints = FALSE,
      all.bin = TRUE,
      use.rw=TRUE)

    # Handle results
    if (prob$status==0){
      clsel = sum(prob$solution)
      auxresultados <- data.frame(auxeblups[,c(ref,traits)], prob$solution)
      resultados <- subset(auxresultados, prob$solution==1)

      if (indsol==0){
        indlinha <- clsel
        clonesout <- resultados[,1]
        colnomes <- c(clsel)
        if (length(resultados) > 3){
          ganhos <- colMeans(resultados[,c(-1, -length(resultados))]*100)
        }else{
          ganhos <- mean(resultados[,c(-1, -length(resultados))]*100)
        }
        gainout <- data.frame(t(ganhos), clsel)
        colnames(gainout) <- c(colnames(resultados[, c(-1, -length(resultados)), drop = FALSE]), "Group.Size")
      }else{
        clonesout <- c(clonesout, resultados[,1], rep(" ", indlinha-clsel))
        colnomes <- c(colnomes, clsel)
        if (length(resultados) > 3){
          ganhos <- colMeans(resultados[,c(-1, -length(resultados))]*100)
        }else{
          ganhos <- mean(resultados[,c(-1, -length(resultados))]*100)
        }
        gainout <- rbind(gainout, c(t(ganhos), clsel))
      }
      indsol <- indsol + 1
    }
  }

  # impossible
  if (indsol == 0) {
    message("No possible solution!")
    return(invisible(NULL))
  }

  clonesout <- matrix(clonesout, nrow = indlinha, ncol = indsol)
  clonesout <- as.data.frame(clonesout)
  colnames(clonesout) <- colnomes


  auxsum <- data.frame(auxsum, row.names = 1)
  auxsum$MaxGain <- NA_real_
  auxsum$MaxGroup <- NA_integer_
  common_names <- intersect(colnames(gainout), rownames(auxsum))

  for (name in common_names) {
    max_val <- max(gainout[[name]], na.rm = TRUE)
    auxsum[name, "MaxGain"] <- max_val
    rows_with_max <- which(gainout[[name]] == max_val)
    last_row <- utils::tail(rows_with_max, 1)
    auxsum[name, "MaxGroup"] <- gainout$Group.Size[last_row]
  }

  if (ncol(clonesout) != (clmax - clmin +1)) {
    cat("No possible solution was found for some of the requested group sizes.")
  }
  result <- list(
    gain = gainout,
    selected = clonesout,
    n_traits = length(traits),
    n_constraints = length(const),
    overview = auxsum
  )
  class(result) <- "polyresult"
  return(result)

}

