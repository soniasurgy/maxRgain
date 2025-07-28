#' Maximum possible gain
#'
#' The maximum possible is the mean of the EBLUPs of the genotypic effects of the best n clones in a given trait, as a percentage of the overall mean. This function calculates the maximum possible gain achieved in the specified traits.
#'
#' @inheritParams polyclonal
#' @note The order of the traits must be consistent across `traits`, `meanvec`, and `criteria`.
#' @returns
#' A list with the following components:
#'  -   `gain`  with the gains of the several traits in each dimension
#'  -   `selected`  with the reference of the clones selected in the group of each dimension
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). \doi{10.1007/s00122-025-04885-0}
#' @export
#' @examples
#' mymeanvec <- c(pa = 12.760)
#' mytraits <- c("pa")
#' maxpos <- rmaxp(
#'    traits = mytraits,
#'    clmin = 7,
#'    clmax = 20,
#'    meanvec = mymeanvec,
#'    data = Gouveio
#'    )
#' maxpos
rmaxp <- function(traits, ref = NULL, clmin = 2, clmax,  meanvec = NULL, criteria = NULL, data)
{
  first <- TRUE
  for (trt in traits){
    mvc <- meanvec[trt]
    crt <- criteria[trt]
    maxpoly <- polyclonal(traits = trt, ref = ref, clmin = clmin, clmax = clmax, meanvec = mvc, criteria = crt, data = data)
    maxp <- as.data.frame(maxpoly$selected)
    gainp <- as.data.frame(maxpoly$gain)
    maxp_lastcol <- length(maxp)
    newclone <- c(maxp[1:clmin, maxp_lastcol])
    for (j in (maxp_lastcol - 1):1){
      newclone <- c(newclone, setdiff(maxp[,j], maxp[,j+1]))
    }
    if (first){
      first <- FALSE
      listclone <- data.frame(temp = newclone)
      listgain <- data.frame(temp = gainp[,2])
      names(listclone)[names(listclone) == "temp"] <- trt
      names(listgain)[names(listgain) == "temp"] <- trt
    }else{
      listclone <- data.frame(listclone, temp = newclone)
      listgain <- data.frame(listgain, temp = gainp[,2])
      names(listclone)[names(listclone) == "temp"] <- trt
      names(listgain)[names(listgain) == "temp"] <- trt
    }
  }
  listclone$Entry <- c(rep(as.character(clmin), clmin), rep("",(clmax-clmin)))
  for (i in (clmin+1):clmax) {
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
#' The maximum admissible genetic gain in one trait that can be achieved without decreasing any of the other traits, as a percentage of the overall mean. This function calculates the maximum admissible gains achieved in the specified traits.
#'
#' @inheritParams polyclonal
#' @param  constraints Vector with traits to which constraints apply. If omitted, all except `ref` are used.
#' @note
#' The order of traits must be consistent across `traits`, `constraints`, `meanvec`, and `criteria`.
#' Both `meanvec` and `criteria` must include values for all traits specified in `traits` and `constraints`.
#' If `constraints` is omitted, all traits in the dataset are considered; in that case, `meanvec` and `criteria` must provide values for all of them.
#' @returns
#' A list with the following components:
#'  -   `gain`  with the gains of the several traits in each dimension
#'  -   `selected_<trait>`  with the reference of the clones selected in the group of each dimension in each trait
#' @references Surgy, S., Cadima, J. & Gonçalves, E. Integer programming as a powerful tool for polyclonal selection in ancient grapevine varieties. Theor Appl Genet 138, 122 (2025). \doi{10.1007/s00122-025-04885-0}
#' @export
#' @examples
#' mymeanvec <- c(yd = 3.517, pa = 12.760, ta = 4.495, ph = 3.927, bw = 1.653)
#' mytraits <- c("yd", "pa")
#' mycriteria <- c(yd = 1, pa = 1, ta = 1, ph = -1, bw = -1)
#' maxadm <- rmaxa(
#'    traits = mytraits,
#'    clmin = 7,
#'    clmax = 20,
#'    meanvec = mymeanvec,
#'    criteria = mycriteria,
#'    data = Gouveio
#'    )
#' maxadm
rmaxa <- function(traits, ref = NULL, clmin = 2, clmax, constraints = NULL, meanvec = NULL, criteria = NULL, data)
{
  selected_list <- list()

  if (!is.null(constraints)){
    auxlength <- length(constraints)

    relmaxa <- c( rep(">=", auxlength))
    rhsmaxa <- c( rep(0, auxlength))
    ctr <- data.frame(a = constraints, b = relmaxa, c = rhsmaxa)
  }else{
    ctr <- all_zero(ref, data)
  }

  first <- TRUE
  for (trt in traits){
    mvc <- meanvec
    maxapoly <- polyclonal(traits = trt, ref = ref, clmin = clmin, clmax = clmax, dmg = ctr, meanvec = mvc, criteria = criteria, data = data)
    maxa <- as.data.frame(maxapoly$selected)
    gaina <- as.data.frame(maxapoly$gain)
    if (first){
      first <- FALSE
      listgain <- data.frame(temp = gaina[,2])
      names(listgain)[names(listgain) == "temp"] <- trt
    }else{
      listgain <- data.frame(listgain, temp = gaina[,2])
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
