###Summary Functions
#' Summarizing Partitioned Diversity
#'
#' @rdname summary.partition
#'
#' @description \code{summary} method for class \code{"partition"}
#'
#' @param object an object of class \code{"partition"}, a result of a call to
#'     \code{\link{partition}}
#' @param x an obect of class \code{"summary_parition"}, a result of a call to
#'     \code{\link{summary.partition}}
#' @param p.value a character; if \code{"one-sided"}, the p-values are interpreted at < observed and significance at p = 0.05.
#'     if \code{"two-sided"}, the p-values are interpreted at both > and < observed and significance at p = 0.025.
#' @param ... additional arguments affecting the summary produced
#'
#' @return The observed and mean expected values of alpha and beta diversity, as
#'     well as the the p-values (number of randomizations < observed diversity)
#'     of each diversity level and interpretations of these. Values are rounded to number of places given by \eqn{1/# Simulations}
#' @return \item{Test}{\code{"ind"} or \code{"sample"} randomization as specified in \code{partition} function}
#' @return \item{P-value}{\code{"one-sided"} or \code{"two-sided"} interpretation of results}
#' @return \item{Randomizations}{number of randomizations run as specified in \code{partition} function}
#' @return \item{q}{Hill number as specified in \code{partition} function}
#' @return \item{Gamma}{observed gammma (regional) diversity}
#' @return \item{Observed}{observed lowest level alpha and additive and multiplicative
#'     beta diversities from the species matrix supplied to \code{\link{partition}}.}
#' @return \item{Expected}{mean lowest level alpha and addtive and multiplicative
#'     beta diversities from randomizations of \code{\link{partition}}.}
#' @return \item{P(<Obs)}{probability of an expected value being lower than the observed value based on number of randomizations}
#' @return Significance interpretations depend on \code{p.value}: if \code{"one-sided"}, the p-values are interpreted at
#'     < observed and significance at p = 0.05. If \code{"two-sided"}, the p-values are interpreted at
#'     > and < observed and significance at p = 0.025; with 0.975 being significant probability of observed being less than expected.
#'
#' @details \code{summary.partition} tries to be smart about formating the
#'     observed and expected diversity values and significance.
#' @details P-values are calculated as \eqn{1 - ((# expected values < observed
#'     values)/# Simulations)}.
#' @details Errors will be thrown if \code{p.value} is not equal to either
#'     \code{"one-sided"} or \code{"two-sided"}.
#'
#' @seealso The diversiting partitioning function \code{\link{partition}}.
#'
#' @importFrom utils write.table
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' ## One-tailed p-values
#' summary(partition.obj)
#'
#' ## Two-tailed p-values
#' summary(partition.obj, p.value = "two-sided")
#' }
#' @export
## S3 method for class "partition"
summary.partition = function(object, p.value = "one-sided",...){
  if(!is(object,"partition")){
    stop("'object' must be the set output of the 'partition' function")
  }
  if(!is.character(p.value)){
    stop("'p.value' must be either 'one-sided' or 'two-sided'")
  }
  Test = q = rands = spaces = Gamma = Hyp.add = Hyp.mult = NULL
  Hyps = Obs = Sigdot = Expt = mean.betas.add = mean.betas.mult = NULL
  AddMult = AlpBeta = partab = part_summ = NULL

  Test <- object$Test
  q <- object$q
  rands <- object$Randomizations
  spaces <- paste("%.",(nchar(1/rands)-2),"f",sep="")
  Gamma <- as.numeric(object$Div[1])
  Hyp.add <- object$Hyp[c(TRUE, FALSE)]
  Hyp.mult <- object$Hyp[c(FALSE, TRUE)]
  Hyps <- round(c(Hyp.add,Hyp.mult[-1]),(nchar(format(1/rands, scientific = FALSE)) - 2))
  #Hyps <- sprintf(spaces,Hyps)
  if(p.value == "one-sided") {
    for(i in 1:length(Hyps)){
      if(Hyps[i] == 0){
        Hyps[i] <- sprintf(spaces,1/rands)
        Sigdot[i] <- "***"
      } else if(Hyps[i] <= 0.001){
        Sigdot[i] <- "***"
      } else if(Hyps[i] <= 0.01){
        Sigdot[i] <- "** "
      } else if(0.05 >= Hyps[i]){
        Sigdot[i] <- "*  "
      } else if(Hyps[i] == 1) {
        Hyps[i] <- sprintf(spaces,1)
        Sigdot[i] <- "   "
      } else {
        Sigdot[i] <- "   "
      }
    }
  } else if(p.value == "two-sided") {
    for(i in 1:length(Hyps)){
      if(Hyps[i] == 0){
        Hyps[i] <- sprintf(spaces,1/rands)
        Sigdot[i] <- "-***"
      } else if(Hyps[i] <= 0.001){
        Sigdot[i] <- "-***"
      } else if(Hyps[i] <= 0.01){
        Sigdot[i] <- "-** "
      } else if(0.025 >= Hyps[i]){
        Sigdot[i] <- "-*  "
      } else if(Hyps[i] == 1) {
        Hyps[i] <- sprintf(spaces,1)
        Sigdot[i] <- "+***"
      } else if(Hyps[i] >= 0.990) {
        Sigdot[i] <- "+** "
      } else if(Hyps[i] >= 0.975) {
        Sigdot[i] <- "+*  "
      }else {
        Sigdot[i] <- "    "
      }
    }
  } else {
    stop("'p.value' must be either 'one-sided' or 'two-sided'")
  }

  Hyps <- sprintf("%9s",sprintf("%.4f",round(as.numeric(Hyps),(nchar(format(1/rands, scientific = FALSE)) - 2))))
  Obs <- sprintf("%8s",sprintf("%.3f",round(as.numeric(object$Div[-1]),3)))
  Expt <- mean(as.numeric(object$Rand.Alpha))
  if(object$Test %in% c("SAMPLE", "sample")){
    for(i in 1:length(object$Rand.Beta.Add)){
      mean.betas.add[i] <- mean(as.numeric(object$Rand.Beta.Add[[i]]))
    }
  } else {
    for(i in 1:nrow(object$Rand.Beta.Add)){
      mean.betas.add[i] <- mean(as.numeric(object$Rand.Beta.Add[i, ]))
    }
  }
  Expt[2:(length(mean.betas.add) + 1)] <- mean.betas.add
  mean.betas.mult <- sapply(object$Rand.Beta.Mult, function(DF){
    mean(as.numeric(DF))
  })

  Expt[(length(mean.betas.mult) + 2):(length(mean.betas.mult) * 2 + 1)] <-
    mean.betas.mult
  Expt <- sprintf("%8s",sprintf("%.3f",round(Expt,3)))
  AddMult <- sprintf("%-14s",c("", "Additive",
                               rep("", times = (length(mean.betas.add) - 1)),
                               "Multiplicative",
                               rep("", times = (length(mean.betas.mult) - 1))))
  AlpBeta <- sprintf("%-6s",c("Alpha",
                              paste("Beta",
                                    rep(c(1:length(mean.betas.add)), times = 2))))
  partab <- data.frame(AddMult,AlpBeta)
  colnames(partab) <- c(sprintf("%14s",""), sprintf("%6s",""))
  partab$Observed <- Obs; partab$Expected <- Expt
  partab$"Pr(< Obs)" <- Hyps; partab$sig <- Sigdot; colnames(partab)[6] <- "   "
  part_summ <- list("Method"  = Test,
                    "P-value" = p.value,
                    "Randomizations" = rands,
                    "q"     = q,
                    "Gamma" = Gamma,
                    "Table" = partab,
                    "Side"  = p.value)
  class(part_summ) <- "summary_partition"
  part_summ
}

#' @rdname summary.partition
#' @export
## S3 method for class "summary_partition"
print.summary_partition = function(x,...){
  cat(paste("Method:",x[[1]]),sep="\n")
  cat(paste("P-value:",x[[2]]),sep="\n")
  cat(paste("Randomizations:",x[[3]]),sep="\n")
  cat(paste("q:",x[[4]]),sep="\n")
  cat(paste("Gamma:",round(x[[5]],3)),sep="\n")
  cat("Diversity Partitioning:",sep="\n")
  utils::write.table(x[[6]], row.names=FALSE, quote=FALSE,sep="  ")
  cat("---",sep="\n")
  if(x[[7]] == "one-sided"){
    cat("Signif. codes: 0.001 '***' 0.01 '**' 0.05 '*'")
  } else {
    cat("Signif. codes: 0.001 '+***' 0.010 '-**' 0.025 '-*'",sep="\n")
    cat("Signif. codes: 1.000 '-***' 0.990 '+**' 0.975 '+*'")
  }
}
