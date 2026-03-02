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
  cat("Maximum possible gains for each trait correspond to independently selected groups; ")
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

  invisible(x)
}

#' Print method for output_rmaxa objects
#'
#' @param x An object of class \code{output_rmaxa}
#' @param ... Additional arguments passed to methods (currently unused)
#' @return Invisibly returns the input object \code{x}
#' @rdname output_rmaxa
#' @export
print.output_rmaxa <- function(x, ...) {
  cat("Maximum admissible gains for each trait correspond to independently selected groups; ")
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

#' Update method for polyclonal results
#'
#' Allows updating any argument of a previous \code{polyclonal()} call without rewriting the full call.
#' Defaults are preserved for arguments not provided.
#'
#' @param object A \code{polyresult} object returned by \code{polyclonal()}.
#' @param ... Arguments to update in the original \code{polyclonal()} call.
#'
#' @return A new \code{polyresult} object with updated arguments applied.
#'
#' @examples
#' #' Original call
#'
#' selections <- polyclonal(traits = c("yd", "pa", "ta", "ph", "bw"),
#'                   clmin = 7, clmax = 15,
#'                   dmg = "base",
#'                   meanvec = c(yd = 3.517, pa = 12.760, ta = 4.495, ph = 3.927, bw = 1.653),
#'                   criteria = c(yd = 1, pa = 1, ta = 1, ph = -1, bw = -1),
#'                   data = Gouveio)
#'
#' selections
#' # Update clmax
#' selupdate1 <- update(selections, clmax = 10)
#'
#' selupdate1
#'
#' # Update dmg
#' selupdate2 <- update(selupdate1, dmg = data.frame(
#'                      lhs = c("yd", "pa", "ta", "ph", "bw"),
#'                      rel = c(">=", ">=", ">=", ">=", ">="),
#'                      rhs = c(25, 0, 0, 0, 0)
#'                      )
#'                    )
#'
#' selupdate2
#'
#' @export
update.polyresult <- function(object, ...) {
  # Recupera a chamada original
  old_call <- object$call

  # Obtem todos os argumentos da função com defaults
  all_args <- as.list(formals(polyclonal))

  # Preenche com os argumentos usados na chamada original
  for (arg in names(old_call)[-1]) {  # ignora o nome da função
    all_args[[arg]] <- old_call[[arg]]
  }

  # Argumentos novos em ...
  dots <- list(...)

  # Se o utilizador forneceu dmg, força a remoção de qualquer valor antigo
  if ("dmg" %in% names(dots)) {
    all_args$dmg <- NULL
  }

  # Agora combina com os novos argumentos
  new_args <- utils::modifyList(all_args, dots)

  # Chama a função original com os novos argumentos
  do.call("polyclonal", new_args)
}
