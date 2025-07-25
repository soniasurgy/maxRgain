#' Print method for polyresult objects
#'
#' Prints the predicted genetic gains and selected genotypes.
#'
#' @param x An object of class \code{polyresult}
#' @param ... Additional arguments passed to methods (currently unused)
#' @return Invisibly returns the input object \code{x}
#' @rdname polyresult
#' @export
print.polyresult <- function(x, ...) {
cat("Predicted genetic gains as a % of the overall mean\n")
cat("\n$gain \n")
gain_df <- x$gain
rownames(gain_df) <- NULL
print(format(gain_df, justify = "left"), row.names = FALSE)
cat("\n")

cat("\nSelected genotypes (per group size)\n")
cat("\n$selected\n")
selected_df <- x$selected
rownames(selected_df) <- NULL
print(format(selected_df, justify = "left"), row.names = FALSE)

invisible(x)

}

#' Summary method for polyresult objects
#'
#' Provides a summary of the results returned by \code{polyclonal()}.
#'
#' @param object An object of class \code{polyresult}
#' @param ... Further arguments passed to methods (currently unused)
#' @return Invisibly returns the input object \code{object}
#' @rdname polyresult
#' @export
summary.polyresult <- function(object, ...) {
  cat("Summary of Selection Results\n")
  cat("-----------------------------------\n")

  # Number of groups
  if (!is.null(object$selected)) {
    cat("\nNumber of groups selected:", ncol(object$selected), "\n")
    cat("Group sizes:", paste(colnames(object$selected), collapse = ", "), "\n")
  } else {
    cat("\nNo selected individuals.\n")
  }

  # Mostrar restrições, se existirem
  if (!is.null(object$overview)) {
    cat("\nOverview\n")
    print(object$overview)
  }

  invisible(object)
}

#' Print method for output_rmaxp objects
#'
#' @param x An object of class \code{output_rmaxp}
#' @param ... Additional arguments passed to methods (currently unused)
#' @return Invisibly returns the input object \code{x}
#' @rdname output_rmaxp
#' @export
print.output_rmaxp <- function(x, ...) {
  cat("Maximum possible gains for each trait correspond to independently selected groups;\n")
  cat("thus, gains are not directly comparable and should be interpreted separately.\n")
  cat("See 'selected' for group details\n\n")

  cat("Predicted genetic gains as a % of the overall mean\n\n")
  cat("$gain\n")
  gain_df <- x$gain

  rownames(gain_df) <- NULL
  print(format(gain_df, justify = "left"), row.names = FALSE)
  cat("\n")

  cat("\nSelected genotypes (per group)\n")
  cat("\n$selected\n")
  selected_df <- x$selected
  rownames(selected_df) <- NULL
  print(format(selected_df, justify = "left"), row.names = FALSE)

}

' Print method for output_rmaxa objects
#'
#' @param x An object of class \code{output_rmaxa}
#' @param ... Additional arguments passed to methods (currently unused)
#' @return Invisibly returns the input object \code{x}
#' @rdname output_rmaxa
#' @export
print.output_rmaxa <- function(x, ...) {
  cat("Maximum admissible gains for each trait correspond to independently selected groups;\n")
  cat("thus, gains are not directly comparable and should be interpreted separately.\n")
  cat("See 'selected' for group details\n\n")

  cat("Predicted genetic gains as a % of the overall mean\n")
  cat("\n$gain \n")
  gain_df <- x$gain

  rownames(gain_df) <- NULL
  print(format(gain_df, justify = "left"), row.names = FALSE)
  cat("\n")

  for (name in names(x)) {
    if (startsWith(name, "selected_")) {
      cat("\n$",name,"\n")
      print(x[[name]], row.names = FALSE)
    }
  }
  invisible(x)
}
