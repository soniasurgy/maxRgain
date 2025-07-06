#' Maximum possible gain
#'
#' The maximum possible is the mean of the EBLUPs of the genotypic effects of the best n clones in a given trait, as a percentage of the overall- This function calculates the maximum possible gain achieved in the specified trait for groups from 7 to 20 clones
#'
#' @inheritParams polyclonal
#' @returns
#' A list with the following components:
#'  -   `gain`  with the gains of the several traits in each dimension
#'  -   `selected`  with the reference of the clones selected in the group of each dimension
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). https://doi.org/10.1007/s00122-025-04885-0
#' @export
#' @examples
#' mymeanvec <- c(pa = 12.760)
#' mytraits <- c("pa")
#' maxpos <- Rmaxp(
#'    traits = mytraits,
#'    meanvec = mymeanvec,
#'    data = Gouveio
#'    )
#' maxpos
Rmaxp <- function(traits, ref = NULL, meanvec = NULL, criteria = NULL, data)
{
  for (i in 1:length(traits)){
    trt <- traits[i]
    mvc <- meanvec[i]
    maxpoly <- polyclonal(traits = trt, ref = ref, clmin = 7, clmax = 20, meanvec = mvc, criteria = criteria, data = data)
    maxp <- as.data.frame(maxpoly$selected)
    gainp <- as.data.frame(maxpoly$gain)
    newclone <- c(maxp[1:7, 14])
    for (j in 13:1){
      newclone <- c(newclone, setdiff(maxp[,j], maxp[,j+1]))
    }
    if (i == 1){
      listclone <- data.frame(temp = newclone)
      listgain <- data.frame(temp = gainp[,1])
      names(listclone)[names(listclone) == "temp"] <- trt
      names(listgain)[names(listgain) == "temp"] <- trt
    }else{
      listclone <- data.frame(listclone, temp = newclone)
      listgain <- data.frame(listgain, temp = gainp[,1])
      names(listclone)[names(listclone) == "temp"] <- trt
      names(listgain)[names(listgain) == "temp"] <- trt
    }
  }
  listclone$Entry <- c(rep("7", 7), rep("",13))
  for (i in 8:20) {
    listclone$Entry[i] <- as.character(i)
  }
  listgain <- data.frame(Group.Size = gainp$Group.Size, listgain)
  result <- list(
    gain = listgain,
    selected = listclone
  )
  class(result) <- "output_rmaxp"
  return(result)
}


#' Maximum admissible gain
#'
#' The maximum genetic gain in one trait that can be achieved without decreasing any of the other traits, as a percentage of the overall mean. This function calculates the maximum admissible gains achieved in the specified trait for polyclonal groups from 7 to 20 clones.
#'
#' @inheritParams polyclonal
#' @param  constraints A vector with the traits for the constraints. If not provided, all the traits in data are assumed.
#' @returns
#' A list with the following components:
#'  -   `gain`  with the gains of the several traits in each dimension
#'  -   `selected_<trait>`  with the reference of the clones selected in the group of each dimension in each trait
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). https://doi.org/10.1007/s00122-025-04885-0
#' @export
#' @examples
#' mymeanvec <- c(yd = 3.517, pa = 12.760, ta = 4.495, ph = 3.927, bw = 1.653)
#' mytraits <- c("yd", "pa")
#' mycriteria <- c(yd = 1, pa = 1, ta = 1, ph = -1, bw = -1)
#' maxadm <- Rmaxa(
#'    traits = mytraits,
#'    meanvec = mymeanvec,
#'    criteria = mycriteria,
#'    data = Gouveio
#'    )
#' maxadm
Rmaxa <- function(traits, ref = NULL, constraints = NULL, meanvec = NULL, criteria = NULL, data)
{
  selected_list <- list()
  if (is.null(constraints)){
    auxlength <- length(data[,-1])
    if (is.null(ref)){
      auxconst <- names(data[,-1])
    }else{
      numcol <- which(names(data) == ref)
      auxconst <- names(data[,-numcol])
    }
  }else{
    auxlength <- length(constraints)
    auxcont <- constraints
  }
  relmaxa <- c( rep(">=", auxlength))
  rhsmaxa <- c( rep(0, auxlength))
  ctr <- data.frame(a = auxconst, b = relmaxa, c = rhsmaxa)
  for (i in 1:length(traits)){
    trt <- traits[i]
    mvc <- meanvec
    maxapoly <- polyclonal(traits = trt, ref = ref, clmin = 7, clmax = 20, dmg = ctr, meanvec = mvc, criteria = criteria, data = data)
    maxa <- as.data.frame(maxapoly$selected)
    gaina <- as.data.frame(maxapoly$gain)
    if (i == 1){
      listgain <- data.frame(temp = gaina[,1])
      names(listgain)[names(listgain) == "temp"] <- trt
    }else{
      listgain <- data.frame(listgain, temp = gaina[,1])
      names(listgain)[names(listgain) == "temp"] <- trt
    }
    name_obj <- paste0("selected_", trt)
    selected_list[[name_obj]] <- maxa
  }
  listgain <- data.frame(Group.Size = gaina$Group.Size, listgain )
  result<- c(list(
    gain = listgain),
    selected_list
  )
  class(result) <- "output_rmaxa"
  return(result)
}


#' Base Situation
#'
#'Situation where the objective function defined by all traits is maximized and the desired minimum gain for each trait is 0 (all the constraints >= 0).
#'
#' @inheritParams polyclonal
#' @param constraints A vector with the traits in the constraints.
#' @returns
#' A list with the following components:
#'  -   `gain`  with the gains of the several traits in each dimension
#'  -   `selected`  with the reference of the clones selected in the group of each dimension
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). https://doi.org/10.1007/s00122-025-04885-0
#' @export
#' @examples
#' mymeanvec <- c(yd = 3.517, pa = 12.760, ta = 4.495, ph = 3.927, bw = 1.653)
#' mytraits <- c("yd", "pa", "ta", "ph", "bw")
#' myconst <- c("yd", "pa", "ta", "ph", "bw")
#' mycriteria <- c(yd = 1, pa = 1, ta = 1, ph = -1, bw = -1)
#' bassit <- BaseSituation(
#'    traits = mytraits,
#'    constraints = myconst,
#'    meanvec = mymeanvec,
#'    criteria = mycriteria,
#'    data = Gouveio
#'    )
#' bassit
BaseSituation <- function(traits, ref = NULL, constraints = NULL, meanvec = NULL, criteria = NULL, data)
{
  if (is.null(constraints)){
    auxlength <- length(data[,-1])
    if (is.null(ref)){
      auxconst <- names(data[,-1])
    }else{
      numcol <- which(names(data) == ref)
      auxconst <- names(data[,-numcol])
    }
  }else{
    auxlength <- length(constraints)
    auxconst <- constraints
  }
  relbase <- c( rep(">=", auxlength))
  rhsbase <- c( rep(0, auxlength))
  ctr <- data.frame(a = auxconst, b = relbase, c = rhsbase)
  return(polyclonal(traits = traits, ref = ref, clmin = 7, clmax = 20, dmg = ctr, meanvec = meanvec, criteria = criteria, data = data))
}
