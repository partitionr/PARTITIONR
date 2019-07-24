## na.q
## function for getting past issue of 0^0 = 1; therefore, need to use
na.q <- function(x,w){
  if(is.na(x)){
    return(0)
  } else if(x == 0){
    return(0)
  }
  else x ^ w
}

## indipart
## function for individual randomization partitioning

indipart <- function(species.num, h.level, n1, sim.rand, q, sp.prop){

  spp.vec <- c()
  rvec <- rand.betas.add <- rand.betas.mult <-  rand.alpha <- pval.add <- NULL
  al.ob <- pval.mult <- hyp.alpha <- al.ob2 <- hyp.gamma <- NULL

  for(i in 1:(length(species.num)))  {
    spp.vec <- c(spp.vec, rep(colnames(species.num[i]),
                              each = sum(species.num[ ,i])))
  }

  if(n1 == 1){
    hyp.alpha <- array(0, sim.rand)
    hyp.gamma <- array(0, sim.rand)
    rand.betas.add <- array(0, sim.rand)
    rand.betas.mult <- array(0, sim.rand)

    if(q != 1){
      for (i in 1:sim.rand) {                                          # For each krand simulations, do the following:
        rvec <- split(spp.vec, sample(length(rowSums(sp.prop[[1]])),
                                      length(spp.vec), replace = TRUE)) # Randomize INDIVIDUALs into Sites
        al.ob <- plyr::rbind.fill(
          lapply(
            lapply(rvec, as.data.frame),
            function(DF1) {
              data.frame(
                prop.table(
                  as.matrix(as.data.frame.matrix(t(table(DF1)))), 1))
            })
        )

        al.ob2 <- data.frame(apply(al.ob, c(1,2), function(DF3){na.q(DF3,w=q)}))

        hyp.alpha[i] <- ((sum(rowSums(al.ob2, na.rm = TRUE)*
                                (1 / length(rownames(al.ob2)))))^(1 / (1 - q)))
        hyp.gamma[i] <- (sum(((1 / length(rownames(al.ob)))*
                                colSums(al.ob, na.rm = TRUE))^q))^(1 / (1 - q))
        rand.betas.add[i] <- hyp.gamma[i] - hyp.alpha[i]
        rand.betas.mult[i] <- hyp.gamma[i]/hyp.alpha[i]
      }
      rand.alpha <- hyp.alpha[i]
    }

    if(q == 1){
      for (i in 1:sim.rand) {                                          # For each krand simulations, do the following:
        rvec <- split(spp.vec, sample(length(rowSums(sp.prop[[1]])),
                                      length(spp.vec), replace = TRUE)) # Randomize INDIVIDUALs into Sites
        al.ob <- plyr::rbind.fill(
          lapply(
            lapply(
              lapply(rvec, as.data.frame),
              function(DF1) {
                data.frame(
                  prop.table(
                    as.matrix(as.data.frame.matrix(t(table(DF1)))), 1))
              }), function(DF2) {
                data.frame(apply(DF2, c(1,2), function(DF3){na.q(DF3,w=q)}))
              }
          )
        )
        al.ob2 <- al.ob*log(al.ob)
        hyp.alpha[i] <- exp(sum(rowSums(al.ob2, na.rm = TRUE)*
                                  (-1 / length(rowSums(al.ob2)))))
        hyp.gamma[i] <- exp(-sum((1 / length(rowSums(al.ob)))*
                                   colSums(al.ob, na.rm = TRUE)*
                                   log((1 / length(rowSums(al.ob)))*
                                         colSums(al.ob, na.rm = TRUE))))

        rand.betas.add[i] <- hyp.gamma[i] - hyp.alpha[i]
        rand.betas.mult[i] <- hyp.gamma[i]/hyp.alpha[i]
      }
      rand.alpha <- hyp.alpha
    }

    Output <- list("Div"            = 0,
                   "Hyp"            = 0,
                   "Rand.Beta.Add"  = rand.betas.add,
                   "Rand.Beta.Mult" = rand.betas.mult,
                   "Rand.Alpha"     = rand.alpha)
    return(Output)     # when calling function, return all three indices
  }

  else if(n1 != 1){
    for(i in 1:length(sp.prop)){
      hyp.alpha[[i]] = array(0,length(rowSums(sp.prop[[i]])))
      rand.betas.mult[[i]]=array(0, sim.rand)
    }
    rand.betas.add <- data.frame(matrix(NA, nrow = length(sp.prop), ncol = sim.rand))

    if(q != 1){
      rvec <- replicate(sim.rand,list())
      for(i in 1:sim.rand){
        for(h in 1:length(sp.prop)){
          rvec[[i]][[h]] <- split(spp.vec, sample(length(rowSums(sp.prop[[h]])),
                                                  length(spp.vec), replace = TRUE)) # Randomize INDIVIDUALs into Sites
        }
      }

      al.ob <-  rvec %>% purrr::modify_depth(1, function(DF){
        lapply(
          lapply(
            lapply(DF, function(sublist) {
              lapply(sublist, as.data.frame)
            }), function(sublist) {
              lapply(sublist, function(DF1){
                data.frame(
                  prop.table(
                    as.matrix(as.data.frame.matrix(t(table(DF1)))), 1)
                )
              })
            }
          ), plyr::rbind.fill)
      })

      al.ob2 <- al.ob %>% purrr::modify_depth(2, function(DF){
        data.frame(apply(DF, c(1,2), function(DF3){na.q(DF3,w=q)}))
      })

      hyp.alpha <- al.ob2 %>% purrr::modify_depth(1,function(DF1){
        sapply(DF1,function(DF2){(sum(rowSums(DF2, na.rm = TRUE)*
                                        (1 / length(rownames(DF2)))))^
            (1/(1 - q))})
      })

      hyp.gamma <- al.ob %>% purrr::modify_depth(1,function(DF1){
        (sum(((1 / length(rowSums(DF1[[1]]))) *
                colSums(DF1[[1]], na.rm = TRUE)) ^ q)) ^ (1 / (1 - q))
      })

      for(i in 1:sim.rand){
        for(h in 1:(length(sp.prop) - 1)){
          rand.betas.add[h,i] <- hyp.alpha[[i]][(h + 1)] - hyp.alpha[[i]][h]
        }
        rand.betas.add[length(sp.prop),i] <- hyp.gamma[[i]] - hyp.alpha[[i]][1] -
          sum(rand.betas.add[c(1:(length(sp.prop) - 1)), i])
        for(h in 1:(length(sp.prop) - 1)){
          rand.betas.mult[[h]][i] <- hyp.alpha[[i]][(h + 1)] / hyp.alpha[[i]][h]
        }
        rand.betas.mult[[(length(sp.prop))]][i] <- hyp.gamma[[i]]/hyp.alpha[[i]][(length(sp.prop))]

        rand.alpha[i] <- hyp.alpha[[i]][1]

      }
    }

    if(q == 1){
      rvec <- replicate(sim.rand,list())
      for(i in 1:sim.rand){
        for(h in 1:length(sp.prop)){
          rvec[[i]][[h]] <- split(spp.vec, sample(length(rowSums(sp.prop[[h]])),
                                                  length(spp.vec), replace = TRUE)) # Randomize INDIVIDUALs into Sites
        }
      }
      al.ob <-  rvec %>% purrr::modify_depth(1, function(DF){
        lapply(
          lapply(
            lapply(DF, function(sublist) {
              lapply(sublist, as.data.frame)
            }), function(sublist) {
              lapply(sublist, function(DF1){
                data.frame(
                  prop.table(
                    as.matrix(as.data.frame.matrix(t(table(DF1)))), 1)
                )
              })
            }
          ), plyr::rbind.fill)
      })

      al.ob2 <- al.ob %>% purrr::modify_depth(2, function(DF){
        DF * log(DF)
      })

      hyp.alpha <- al.ob2 %>% purrr::modify_depth(1,function(DF1){
        sapply(DF1,function(DF2){exp(sum(rowSums(DF2, na.rm = TRUE)*
                                           (-1/length(rownames(DF2)))))
        })
      })

      hyp.gamma <- al.ob %>% purrr::modify_depth(1,function(DF1){
        exp(-sum((1/length(rownames(DF1[[length(sp.prop)]])))*
                   colSums(DF1[[length(sp.prop)]], na.rm = TRUE)*
                   log((1/length(rownames(DF1[[length(sp.prop)]])))*
                         colSums(DF1[[length(sp.prop)]], na.rm = TRUE))))
      })

      for(i in 1:sim.rand){
        for(h in 1:(length(sp.prop) - 1)){
          rand.betas.add[h,i] <- hyp.alpha[[i]][(h + 1)] - hyp.alpha[[i]][h]
        }
        rand.betas.add[length(sp.prop),i] <- hyp.gamma[[i]] - hyp.alpha[[i]][1] -
          sum(rand.betas.add[c(1:(length(sp.prop) - 1)), i])
        for(h in 1:(length(sp.prop) - 1)){
          rand.betas.mult[[h]][i] <- hyp.alpha[[i]][(h + 1)] / hyp.alpha[[i]][h]
        }
        rand.betas.mult[[(length(sp.prop))]][i] <- hyp.gamma[[i]]/hyp.alpha[[i]][(length(sp.prop))]

        rand.alpha[i] <- hyp.alpha[[i]][1]

      }
    }

    Output <- list("Div"            = 0,
                   "Hyp"            = 0,
                   "Rand.Beta.Add"  = rand.betas.add,
                   "Rand.Beta.Mult" = rand.betas.mult,
                   "Rand.Alpha"     = rand.alpha)
    return(Output)     # when calling function, return all three indices
  }

}


## samplpart
## function for sample randomization partitioning

samplpart <- function(species.num, h.level, n1, sim.rand, q, factors, sp.use){

  spdat <- data.frame(species.num)

  rand.beta.add.samp <- rand.beta.mult.samp <- f.rand.alpha <- NULL
  real.rand.betas.add.samp <- real.rand.betas.mult.samp <- NULL

  for(i in 2:(n1)){
    if((i+1) == n1){
      prene <- as.factor(factors[,h.level[n1]])
    } else if(i == n1){
      prene <- as.factor(factors[,h.level[n1]])
    } else {
      prene <- as.factor(apply(factors[,c((i+1):n1)], 1, paste, collapse="."))
    }

    nummer <- dplyr::count_(as.data.frame(prene), "prene")
    annoyance <- factors[ , c((i-1):n1)]

    y <- list(); yy <- list(); yyg <- list(); yy.summ <- NULL; num <- list()
    alpha.vec <- list(); real.rand.alpha <- NULL; rand.alpha <- list()
    f.rand.alpha[[i]] <- array(0, 1)
    real.rand.betas.add.samp[[i]] <- array(0, 1)
    real.rand.betas.mult.samp[[i]] <- array(0, 1)

    if(i != n1){
      for(f in 1:length(nummer$n)) {
        gravy <-NULL
        gravy <- prene == nummer$prene[f]
        speck <- NULL
        sprite <- data.frame(matrix(nrow = 0, ncol = ncol(spdat)))
        for(g in 1:length(gravy)) {
          if(gravy[g] == TRUE) {
            sprite[g, ] <- spdat[g, ]
          } else {
            sprite[g, ] <- rep(NA, ncol(spdat))
          }
        }

        speck <- c(speck, sprite)

        speck <- as.data.frame(matrix(do.call(rbind, speck),
                                      nrow = length(gravy),
                                      ncol = ncol(spdat),
                                      byrow = TRUE))
        speck <- speck[stats::complete.cases(speck), ]
        rows <- rownames(speck)
        t1 <- sp.use[rows, ]                                        # Finally! This is the hierarchical level in which everything remains the same
        t1 <- droplevels(t1)

        fixt <- sapply(t1, is.numeric)

        spdat_new <- data.frame(t1[ , fixt])

        facts <- factors[rows, ]

        if((i - 1) != 1){
          factim1 <- interaction(c(facts[ , c((i-1):n1)]))
          spdat_new$im1level <- factim1
        }

        facts <- droplevels(facts[ , which(names(facts) %in%
                                             colnames(annoyance[ ,c(2:length(annoyance))]))])
        if((i+1) != n1){
          factsip1 <- droplevels(facts[ , which(names(facts) %in%
                                                  colnames(annoyance[ , c(3:length(annoyance))]))])
        } else {
          factsip1 <- droplevels(data.frame(facts[ , 2]))
        }

        y <- split(spdat_new, factsip1)               # and the level in which they will be placed (i+1)
        yy[[f]] <- y

        if((i - 1) != 1){
          yy.summ[[f]] <- data.frame(yy[[f]])
          colnames(yy.summ[[f]])[ncol(yy.summ[[f]])] <- "im1level"
          yyg[[f]] <- data.frame(yy.summ[[f]] %>% dplyr::group_by_(~im1level) %>%
                                   dplyr::summarise_all(dplyr::funs(sum)))
          yyg[[f]] <- yyg[[f]][ , -1]
        } else {
          yyg[[f]] <- yy[[f]]
        }
      }

      num.list <- rand.spdat <- rand.spdat.i <- rand.sp.grp.i <- NULL
      rand.prop.ip1 <- rand.alpha.ip1 <- rand.prop.i <- rand.spdat.ip1 <- NULL
      rand.sp.grp.ip1 <- rand.prop.i.2 <- rand.prop.ip1.2 <- rand.alpha.i <- NULL
      rand.alphas <- rand.beta.add <- rand.beta.mult <- rand.prop <- NULL
      rand.alphas.p1 <- NULL

      for(j in 1:sim.rand){
        num.list <-  lapply(
          lapply(yyg,function(DF1){
            dplyr::sample_n(data.frame(DF1), size = nrow(data.frame(DF1)),
                            replace = FALSE)
          }), function(DF2){
            names(DF2) <- c(1:length(colnames(DF2)))
            DF2
          })
        rand.spdat <- data.frame(do.call(rbind,num.list))
        rand.spdat.i <- cbind(rand.spdat,
                              interaction(dplyr::distinct_(annoyance)[ , colnames(facts)]))
        colnames(rand.spdat.i)[length(rand.spdat.i)] <- c("ilevel")
        rand.sp.grp.i <- data.frame(rand.spdat.i %>%
                                      dplyr::group_by_(~ilevel) %>%
                                      dplyr::summarise_all(dplyr::funs(sum)))
        rand.prop.i <- prop.table(as.matrix(rand.sp.grp.i[ ,
                                                           c(2:ncol(rand.sp.grp.i))]), 1)
        rand.spdat.ip1 <- suppressWarnings(cbind(rand.spdat,
                                                 interaction(dplyr::distinct_(annoyance)[ , colnames(facts)[2:length(facts)]])))
        colnames(rand.spdat.ip1)[length(rand.spdat.ip1)] <- c("ip1level")
        rand.sp.grp.ip1 <- data.frame(rand.spdat.ip1 %>%
                                        dplyr::group_by_(~ip1level) %>%
                                        dplyr::summarise_all(dplyr::funs(sum)))
        rand.prop.ip1 <- prop.table(as.matrix(rand.sp.grp.ip1[ , c(2:ncol(rand.sp.grp.ip1))]), 1)

        if(q != 1){
          rand.prop.i.2 <- rand.prop.i
          rand.prop.i.2 <- apply(rand.prop.i, c(1, 2), function(DF3){na.q(DF3,w=q)})
          rand.alpha.i <- ((sum(rowSums(rand.prop.i.2, na.rm = TRUE)*
                                  (1/length(rowSums(rand.prop.i.2)))))^(1/(1-q)))

          rand.prop.ip1.2 <- rand.prop.ip1
          rand.prop.ip1.2 <- apply(rand.prop.ip1, c(1, 2), function(DF3){na.q(DF3,w=q)})
          rand.alpha.ip1 <- ((sum(rowSums(rand.prop.ip1.2, na.rm = TRUE)*
                                    (1/length(rowSums(rand.prop.ip1.2)))))^(1/(1-q)))

        } else {
          rand.prop.i.2 <- rand.prop.i
          rand.prop.i.2 <- rand.prop.i*log(rand.prop.i)
          rand.alpha.i <- exp(sum(rowSums(rand.prop.i.2, na.rm = TRUE)*
                                    (-1/nrow(rand.prop.i.2))))
          rand.prop.ip1.2 <- rand.prop.ip1
          rand.prop.ip1.2 <- rand.prop.ip1*log(rand.prop.ip1)
          rand.alpha.ip1 <- exp(sum(rowSums(rand.prop.ip1.2, na.rm = TRUE)*
                                      (-1/nrow(rand.prop.ip1.2))))
        }
        rand.alphas[j] <- rand.alpha.i
        rand.alphas.p1[j] <- rand.alpha.ip1
        rand.beta.add[j] <- rand.alpha.ip1 - rand.alpha.i
        rand.beta.mult[j] <- rand.alpha.ip1 / rand.alpha.i
      }
    } else {
      n1dat <- cbind(spdat, interaction(factors[ , c((n1 - 1), n1)]))
      colnames(n1dat)[ncol(n1dat)] <- "ilevel"
      yyg <- data.frame(n1dat %>%
                          dplyr::group_by_(~ilevel) %>%
                          dplyr::summarise_all(dplyr::funs(sum)))
      yyg <- yyg[ , -1]

      num.list <- rand.spdat <- rand.sp.grp <- rand.prop <- rand.prop.2 <- NULL
      rand.alpha <- rand.gamma <- rand.alphas <- rand.beta.add <- NULL
      rand.beta.mult <- NULL
      rand.gammas <- NULL
      for(j in 1:sim.rand){
        num.list = data.frame(0);rand.prop = data.frame(0)
        num.list <- sample_n(data.frame(yyg), size = nrow(data.frame(yyg)),
                             replace = FALSE)
        colnames(num.list) <- c(1:ncol(num.list))
        rand.spdat <- data.frame(num.list)
        rand.spdat <- suppressWarnings(cbind(rand.spdat,
                                             dplyr::distinct_(annoyance)))
        colnames(rand.spdat)[length(rand.spdat)] <- c("ilevel")
        rand.spdat <- rand.spdat[ , -(ncol(rand.spdat) - 1)]
        rand.sp.grp <- data.frame(rand.spdat %>% dplyr::group_by_(~ilevel) %>%
                                    dplyr::summarise_all(dplyr::funs(sum)))
        rand.prop <- prop.table(as.matrix(rand.sp.grp[ ,
                                                       c(2:ncol(rand.sp.grp))]) , 1)

        if(q != 1){
          rand.prop.2 <- rand.prop
          rand.prop.2 <- apply(rand.prop,c(1, 2), function(DF3){na.q(DF3,w=q)})
          rand.alpha <- ((sum(rowSums(rand.prop.2, na.rm = TRUE)*
                                (1/length(rowSums(rand.prop.2)))))^(1/(1 - q)))

          rand.gamma <- (sum(((1 / nrow(rand.prop)) *
                                colSums(rand.prop, na.rm = TRUE)) ^ q)) ^ (1 / (1 - q))


        } else {
          rand.prop.3 <- rand.prop
          rand.prop.3 <- rand.prop*log(rand.prop)

          rand.alpha <- exp(sum(rowSums(rand.prop.3, na.rm = TRUE)*
                                  (-1/nrow(rand.prop.3))))

          rand.gamma <- exp(-sum((1 / nrow(rand.prop)) *
                                   colSums(rand.prop, na.rm=TRUE)*
                                   log((1 / nrow(rand.prop)) *
                                         colSums(rand.prop, na.rm = TRUE))))
        }

        rand.gammas[j] <- rand.gamma
        rand.beta.add[j] <- rand.gamma - rand.alpha
        rand.beta.mult[j] <- rand.gamma / rand.alpha
      }
    }

    if(i == 2){
      rand.alphas1 <- rand.alphas
    }

    rand.beta.add.samp[[i]] <- rand.beta.add
    rand.beta.mult.samp[[i]] <- rand.beta.mult
  }

  Output <- list("Div"            = 0,
                 "Hyp"            = 0,
                 "Rand.Beta.Add"  = rand.beta.add.samp[c(2:n1)],
                 "Rand.Beta.Mult" = rand.beta.mult.samp[c(2:n1)],
                 "Rand.Alpha"     = rand.alphas1)
  return(Output)     # when calling function, return all three indices
}


##partition - defaults of no hypothesis test (returns just calculated Gamma,
##            alpha and beta values), q inde equal to 0 (species richness), and
##            simulations equal to 1000 randomizations
##          - this function calculates both additive AND multiplicative Beta
##            diversity as a default (and cannot be changed)
##          - the first "level" of the sampling heirarchy must be the smallest
##            "sample" unit - in the scope of the Canopy Study, the smallest
##            unit was 12 samples taken from each tree. Ths sampling unit is not
##            ecologically important, and is thus not taken into consideration
##            in the calculation and randomization of the data
##          - described as a "large function" at 806.8 Kb
##          - run times are currently a bit long (4 minutes for individual based
##            and 7 minutes for sample based randomizations of 1152 obs. at 5
##            hierarchical levels with 94 species)
##          - returned objects are of the "partition" class, which include
##            -  $Div : observed diversity measures
##            -  $Hyp : outcomes of the hypothesis testing
##            -  $Rand.Beta.Add
##                    : expected additive Beta diversity
##            -  $Rand.Beta.Mult
##                    : expected multiplicative Beta diversity
##            -  $Rand.Alpha
##                    : expected alpha diversity at the lowest ecologically
##                      important level
##            -  $Test
##                    : "INDIVIDUAL" or "SAMPLE" - set by user, passed to object
##                      to be used in secondary and generic functions
##            -  $q
##                    : numeric - set by user, passed to object
##                      to be used in secondary and generic functions
##            -  $Randomizations
##                    : numeric - set by user, passed to object
##                      to be used in secondary and generic functions

#' Partitioning Diversity
#'
#' @description \code{partition} is used to calculate alpha, beta, and gamma
#'     diversity of both balanced and unbalanced designs. It can be used to
#'     carry out diversity partitioning at both single and multiple scales, as
#'     well as nested designs. \code{partition} can run hypothesis testing based
#'     on distributions generated from two randomizations methods (see below).
#'
#' @param sp Data frame containing a species matrix that includes variables
#'     coding for the different scales by which diversity is to be partitioned.
#'     The columns containing levels by which diversity is to be partitioned
#'     \strong{must} be of class \code{factor}.
#' @param h.level Vector or list containing column names coding for the
#'     different scales by which diversity is to be partitioned. These must be
#'     in increasing scale.
#' @param low.level Integer representing the lowest level of diversity to be
#'     partitioned. This is the lowest ecologically relevant level. Defaults to
#'     \code{1}. Must be \code{>1} when \code{hyp.test = "SAMPLE"}.
#' @param q Integer representing Hill Number q-diversity metrics. \code{0}
#'     represents species richness; \code{1} represents Shannon-diversity; and
#'     \code{2} represents Simpsons-diversity. See Jost (2007) for more
#'     information.
#' @param hyp.test Method of hypothesis testing to be used; for hypothesis
#'     testing: \code{hyp.test = "INDIVIDUAL"} and \code{hyp.test = "SAMPLE"};
#'     for calculation of observed values only: \code{hyp.test = "NONE"}.
#' @param sim.rand Integer representing the number of randomizations to be run
#'     for hypothesis testing; defaults to \code{1000}. If
#'     \code{hyp.test = "NONE"}, \code{sim.rand} can be ignored.
#' @param x object returned by \code{partition} of class \code{"partition"}
#' @param ... further arguments passed to or from other methods
#'
#' @details Diversity partitioning is a method of decomposing a total amount of
#'     diversity (\emph{gamma} - \eqn{\gamma}) into the components of mean
#'     diversity within samples (\emph{alpha} - \eqn{\alpha}) and diversity
#'     among samples (\emph{beta} - \eqn{\beta}). It can be used with a wide
#'     variety of diversity metrics, specifically the Hill Number
#'     \emph{q}-diversity metrics (Jost 2007). \eqn{\gamma} and \eqn{\alpha} are
#'     calculated using equations 3 - 6 of Chao et al. (2012), with equations 5
#'     and 6 used when \code{q = 1} and equations 3 and 4 used in all other
#'     cases.
#' @details The \code{partition} function can be applied to a variety of data
#'     sets, including those with an unbalanced sampling design, substantial
#'     variation in the number of individuals within those samples, and
#'     hierarchical sampling designs (multiple nested levels).
#' @details \code{partition} calculates alpha and beta diversity and uses
#'     randomization to derive expected values of alpha and beta diversity that
#'     would be obtained if individuals (\code{INDIVIDUAL}) or samples
#'     (\code{SAMPLE}) were randomly distributed. This randomization allows for
#'     significance testing of the observed diversity estimates. The statistical
#'      rationale and operational description of \code{INDIVIDUAL}- and
#'     \code{SAMPLE}-based randomization can be found in Crist et al. (2003).
#' @details The \code{partition} function is the R equivalent of the PARTITION
#'     software developed by Crist et al. (2003).
#' @details At the highest sampling level (h), the diversity components are
#'     calculated as follows:
#'     \deqn{Additive: \beta (h) = \gamma - \alpha (h)}
#'     \deqn{Multiplicative: \beta (h) = \gamma / \alpha (h)}
#'     For the lower sampling levels calculated as follows:
#'     \deqn{Additive: \beta (i) = \alpha (i+1) - \alpha (i)}
#'     \deqn{Multiplicative: \beta (i) = \alpha (i+1) / \alpha (i)}
#'
#' @return An object of class \code{"partition"}. This object is a list of data
#'     frames. Calling an object of class \code{"partition"} will print the
#'     \code{$Div} and \code{$Hyp} data frames; these data frames contain the
#'     partitioned data and p-values from significance tests, respectively. For
#'     more information from the object and interpretations of significance
#'     tests, please use the \code{\link{summary.partition}} function.
#' @return \item{\code{$Div}}{observed partitioned diversity}
#' @return \item{\code{$Hyp}}{p-values from significance testing}
#' @return \item{\code{$Rand.Beta.Add}}{expected additive Beta diversity}
#' @return \item{\code{$Rand.Beta.Mult}}{expected multipicative Beta diversity}
#' @return \item{\code{$Rand.Alpha}}{expected alpha diversity at the lowest
#'     ecologically important level (set by \code{low.level})}
#' @return \item{\code{$Test}}{\code{"INDIVIDUAL"}, \code{"SAMPLE"}, or
#'     \code{"NONE"} - set by \code{hyp.test}, passed to object to be used in
#'     support and generic functions}
#' @return \item{\code{$q}}{integer set by \code{q}, passed to object to be used
#'     in support and generic functions}
#' @return \item{\code{$Randomizations}}{integer set by \code{sim.rand}, passed
#'     to object to be used in support and generic functions}
#'
#' @importFrom plyr rbind.fill
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise_all
#' @importFrom dplyr count_
#' @importFrom dplyr sample_n
#' @importFrom dplyr distinct_
#' @importFrom dplyr funs
#' @importFrom purrr modify_depth
#' @importFrom stats complete.cases
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' part.obj <- partition(sp = spiders.spp,
#'                       h.level = c("SAMPLE", "TREESP"),
#'                       low.level = 1,
#'                       q = 0,
#'                       hyp.test = "INDIVIDUAL",
#'                       sim.rand = 1000)
#'
#' print(part.obj)
#' }
#' @export
partition <- function(sp, h.level, low.level = 1, q = 0, hyp.test = "INDIVIDUAL", sim.rand = 1000) {
  if(low.level == 1 & hyp.test == "SAMPLE"){
    stop("Cannot have 1st level be lowest level of analysis with sample-based
         randomizations")
  }
  n1 <- length(h.level)
  s <- NULL; spdat <- NULL; rich <- NULL; num <- NULL; sp.mat <- NULL
  alpha <- NULL; beta.add <- NULL; beta.mult <- NULL;	sp.prop <- NULL
  sp.prop2 <- NULL; sp.prop3 <- NULL
  species.list <- sapply(sp, is.numeric)
  species.num <- sp[ , species.list]
  factors <- data.frame(sp[ , h.level])
  sp.use <- cbind(species.num, factors)
  for(i in low.level:n1){
    spdat[[i]] <- data.frame(species.num)
    if((i) == n1){
      spdat[[i]]$l <- as.factor(factors[ , n1])
    } else {
      spdat[[i]]$l <- as.factor(apply(factors[ , c((i):n1)], 1, paste,
                                      collapse="."))
    }
    rich[[i]] <- data.frame(spdat[[i]] %>% dplyr::group_by_(~l) %>%
                              dplyr::summarise_all(dplyr::funs(sum)))
    num[[i]] <- sapply(rich[[i]], is.numeric)
    sp.mat[[i]] <- rich[[i]][ , num[[i]]]
    sp.prop[[i]] <- prop.table(as.matrix(sp.mat[[i]]), 1)
  }
  if(low.level > 1){
    for(i in rev(1:(low.level - 1))){
      sp.prop[[i]] <- NULL
      spdat[[i]] <- NULL
    }
  }
  if(q != 1){
    alpha <- sapply(
      lapply(sp.prop, function(DF1){
        apply(DF1, c(1, 2),
              function(DF3){na.q(DF3,w=q)})
      }), function(DF2){
        ((sum(rowSums(DF2, na.rm = TRUE) *
                (1/nrow(DF2)))) ^ (1 / (1 - q)))
      })
    gamma <- (sum(((1 / nrow(sp.prop[[length(sp.prop)]])) *
                     colSums(sp.prop[[length(sp.prop)]], na.rm = TRUE)) ^ q)) ^ (1 / (1 - q))

    if(n1 == 1){
      beta.add <- gamma - alpha
      beta.mult <- gamma / alpha
    }
    else if(n1 != 1){
      for(i in 1:(n1 - 1)){
        beta.add[i] <- alpha[(i + 1)] - alpha[i]
      }
      beta.add[length(sp.prop)] <- gamma - alpha[(n1 - 1)]

      for(i in 1:(length(sp.prop) - 1)){
        beta.mult[i] <- alpha[(i + 1)] / alpha[i]
      }
      beta.mult[length(sp.prop)] <- gamma / alpha[(length(alpha))]
    }
  }

  else if(q == 1){
    alpha <- sapply(
      lapply(sp.prop,function(DF1){
        DF1*log(DF1)
      }), function(DF2){
        (exp(sum(rowSums(DF2, na.rm = TRUE) *
                   (-1 / nrow(DF2)))))
      }
    )

    gamma <- exp(-sum((1 / nrow(sp.prop[[length(sp.prop)]])) *
                        colSums(sp.prop[[length(sp.prop)]], na.rm=TRUE)*
                        log((1 / nrow(sp.prop[[length(sp.prop)]])) *
                              colSums(sp.prop[[length(sp.prop)]], na.rm = TRUE))))

    if(n1 == 1){
      beta.add <- gamma - alpha
      beta.mult <- gamma / alpha
    }
    else if(n1 != 1){
      for(i in 1:(length(sp.prop) - 1)){
        beta.add[i] <- alpha[(i + 1)] - alpha[i]
      }
      beta.add[length(sp.prop)] <- gamma - alpha[(n1 - 1)]

      for(i in 1:(length(sp.prop) - 1)){
        beta.mult[i] <- alpha[(i + 1)] / alpha[i]
      }
      beta.mult[length(sp.prop)] <- gamma / alpha[(length(alpha))]
    }
  }
  else {stop("ERROR: q-diversity metric must be a numeric operator")}

  if(hyp.test == "INDIVIDUAL"){
    Output <- indipart(species.num, h.level, n1, sim.rand, q, sp.prop)

    pval.add <- pval.mult <- NULL

    for(i in 1:length(Output$Rand.Beta.Mult)){
      pval.add[i]<- 1-(sum(beta.add[i]>Output$Rand.Beta.Add[i,])/sim.rand)
      pval.mult[i]<- 1-(sum(beta.mult[i]>Output$Rand.Beta.Mult[[i]])/sim.rand)
    }
    pval.alph <- 1-(sum(alpha[1]>Output$Rand.Alpha)/sim.rand)

    hyp.pvals <- matrix(c(pval.alph,pval.add[1:length(sp.prop)], pval.alph,pval.mult[1:length(sp.prop)]),
                        nrow = 2, byrow = TRUE)
    colnames(hyp.pvals) <- c("Alpha 1", paste("Beta", c(1:length(sp.prop))))
    rownames(hyp.pvals) <- c('Additive p-vals', 'Multiplicative p-vals')
    hyp.table <- as.table(hyp.pvals)

    Output$Div            = c("Gamma" = gamma,
                              "Alpha" = alpha[1],
                              "Additive Beta" = beta.add,
                              "Multiplicative Beta" = beta.mult)
    Output$Test           = hyp.test
    Output$q              = q
    Output$Randomizations = sim.rand
    Output$Hyp            = hyp.table

    class(Output) <- c("partition", class(Output))
    return(Output)     # when calling function, return all three indices
  }
  #####################################Sample Hypothesis Test
  if(hyp.test == "SAMPLE"){
    Output <- samplpart(species.num, h.level, n1, sim.rand, q, factors, sp.use)

    pval.alph <- pval.add <- pval.mult <- NULL

    pval.alph <- 1 - (sum(alpha[1]>Output$Rand.Alpha)/sim.rand)

    for(i in 1:length(Output$Rand.Beta.Add)){
      pval.add[i]<- 1 - (sum(beta.add[i]>as.vector(Output$Rand.Beta.Add[[i]]))/(sim.rand))
      pval.mult[i]<- 1 - (sum(beta.mult[i]>as.vector(Output$Rand.Beta.Mult[[i]]))/(sim.rand))
    }
    hyp.pvals <- matrix(c(pval.alph, pval.add[1:length(sp.prop)], pval.alph, pval.mult[1:length(sp.prop)]),
                        nrow = 2, byrow = TRUE)
    colnames(hyp.pvals) <- c("Alpha 1", paste("Beta", c(1:length(sp.prop))))
    rownames(hyp.pvals) <- c('Additive p-vals', 'Multiplicative p-vals')
    hyp.table <- as.table(hyp.pvals)

    Output$Div            = c("Gamma" = gamma,
                              "Alpha" = alpha[1],
                              "Additive Beta" = beta.add,
                              "Multiplicative Beta" = beta.mult)
    Output$Test           = hyp.test
    Output$q              = q
    Output$Randomizations = sim.rand
    Output$Hyp            = hyp.table

    class(Output) <- c("partition", class(Output))
    return(Output)     # when calling function, return all three indices
  }
  #####################################No Hypothesis Test
  if(hyp.test == "NONE"){
    Output <- list("Div" = c("Gamma" = gamma,
                             "Alpha" = alpha[1],
                             "Additive Beta" = beta.add,
                             "Multiplicative Beta" = beta.mult))
  class(Output) <- c("partition", class(Output))
  return(Output)
  }

  else{stop('hyp.test must be set to "NONE", "INDIVIDUAL", or "SAMPLE"')}

  }

#' @rdname partition
#' @export
## S3 method for class "partition"
print.partition = function(x,...)print(x[1])
