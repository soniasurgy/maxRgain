#' @noRd
print.output_ipp <- function(x, ...) {
if (!is.null(x$dmg)) {
  cat("Desired minimum gains (dmg):\n")
  print(unname(x$dmg))
  cat("\n")
}
cat("Predited genetic gains as a % of the overall mean\n")
cat("\n$gain \n")
gain_df <- x$gain
gain_df <- gain_df[, c(ncol(gain_df), 1:(ncol(gain_df)-1))]
rownames(gain_df) <- NULL
print(format(gain_df, justify = "left"), row.names = FALSE)
cat("\n")

cat("\nSelected genotypes (per group size)\n")
cat("\n$selected\n")
selected_df <- x$selected
rownames(selected_df) <- NULL
print(format(selected_df, justify = "left"), row.names = FALSE)
}





#' @noRd
print.output_rmaxp <- function(x, ...) {
  cat("Rmaxp gains for each trait correspond to independently selected groups;\n")
  cat("thus, gains are not directly comparable and should be interpreted separately.\n")
  cat("See 'selected' for group details\n\n")
  
  cat("Predited genetic gains as a % of the overall mean\n\n")
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


#' @noRd
print.output_rmaxa <- function(x, ...) {
  cat("Rmaxa gains for each trait correspond to independently selected groups;\n")
  cat("thus, gains are not directly comparable and should be interpreted separately.\n")
  cat("See 'selected' for group details\n\n")
  
  cat("Predited genetic gains as a % of the overall mean\n")
  cat("\n$gain \n")
  gain_df <- x$gain

  rownames(gain_df) <- NULL
  print(format(gain_df, justify = "left"), row.names = FALSE)
  cat("\n")
  
  for (name in names(x)) {
    if (startsWith(name, "selected_")) {
      cat("\n===>", name, "\n")
      print(x[[name]], row.names = FALSE)
    }
  }
}