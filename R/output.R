output_ipp <- function(indsol, clonesout, indlinha, colnomes, gainout){
  clonesout <- matrix(clonesout, nrow = indlinha, ncol = indsol)
  clonesout <- as.data.frame(clonesout)
  colnames(clonesout) <- c(colnomes)
  return(list(
    gain = gainout,
    selected = clonesout
  ))
}
