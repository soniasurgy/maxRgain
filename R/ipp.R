
ipp <- function(traits, ref, clmin, clmax, const, relation, rhs, data, dmg)
{
  clsel <- 0
  clselmax <- 0
  clselmin <- 0
  indsol <- 0
  objvalue <- 0
  gainout <- NULL
  auxeblups <- data.frame(data)
  #***** Objective function *****
  objdir <- "max"
  auxobj <- auxeblups[, traits, drop = FALSE]
  fobj=0
  for (a in 1:length(auxobj)){
    fobj=fobj+auxobj[a]
  }
  colnames(fobj) <- c("objfun")
  if (is.null(const)){
    const <- c(1)
    constvalue <- data.frame(rep(1,nrow(auxeblups)))
    colnames(constvalue) <- c("var")
  }else{
    auxconst <- auxeblups[const]
    constlength <- length(auxconst)
    constvalue <- data.frame(rep(1,nrow(auxeblups)), auxconst )
    colnames(constvalue) <- c("var", colnames(const))
  }
  if (missing(relation)){
    dirvalue <- c("==")
  }else{
    dirvalue <- c("==",relation)
  }
  for (i in clmax:clmin){
    rhsvalue  <- i
    if (!missing(rhs)){
      for (a in 1:length(rhs)){
        rhsvalue <- c(rhsvalue,  rhs[a]*i/100)
      }
    }
    prob <- lpSolve::lp(objdir, objective.in = fobj, const.mat = constvalue, const.dir=dirvalue, const.rhs=rhsvalue , transpose.constraints = FALSE,  all.bin = TRUE, use.rw=TRUE)
    if (prob$status==0){
      clsel = sum(prob$solution)
      auxresultados <- data.frame(auxeblups[,c(ref,traits)], prob$solution)
      resultados <- subset(auxresultados, prob$solution==1)
      if (indsol==0){
        indlinha <- clsel
      }
      if (indsol==0){
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
  if (indsol==0){
    print("No possible solution!")
  }else{
    clonesout <- matrix(clonesout, nrow = indlinha, ncol = indsol)
    clonesout <- as.data.frame(clonesout)
    colnames(clonesout) <- colnomes

    if (!is.null(dmg)){
      dmg <- data.frame(dmg, row.names = 1)
    }
    result <- list(
      dmg = dmg,
      gain = gainout,
      selected = clonesout
    )
    class(result) <- "output_ipp"
    return(result)
  }
}
