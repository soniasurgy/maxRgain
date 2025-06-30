norm_eblup <- function(dframe, cols, meanvec, criteria){
  eblup_norm <- dframe
  for (col in cols) {
    if (!is.null(meanvec)){
      eblup_norm[[col]] <- eblup_norm[[col]] / meanvec[[col]]
    }
    if (!is.null(criteria)){
      eblup_norm[[col]] <- eblup_norm[[col]] * criteria[[col]]
    }
  }
  return(eblup_norm)
}