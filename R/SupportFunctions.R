#' Plotting Partitioned Diversity
#'
#' @description Graphical display of partitioned data. Graphics include a stacked bar chart
#'     and line plot. Both observed partition and expected partition (i.e., mean of the null
#'     statistical distributions) are shown in both the bar chart and line plot.
#'
#' @param part.obj an object of class \code{"partition"}, a result of a call to
#'     \code{\link{partition}}
#' @param beta.type the method of partition to be plotted; if \code{"add"}, additive diversity will be plotted;
#'     if \code{"mult"}, multiplicative diversity will be plotted.
#' @param plot.type the type of graph; if \code{"bar"}, a stacked bar chart will be plotted.
#'     if \code{"line"}, a line plot will be plotted with 95\% CI error bars.
#'
#' @return A \code{\link{ggplot}} graphic of the observed and expected partition.
#'     \code{"bar"} represents a stacked bar chart with the y-axis representing a
#'     proportion of the gamma (regional) diversity explained by each level. \code{"line"}
#'     represents a line plot of the change in alpha- and beta-diversity of both observed
#'     and expected data. Error bars represent 95\% CIs of null statistical distribution.
#'
#' @details Function \code{partiplot} visualizes the partitioned data from an object of
#'     class \code{partition}.
#' @details The function uses \code{\link[ggplot2]{ggplot}} to visualize the partitioned data.
#'     To save the graphic, use the \code{\link[ggplot2]{ggsave}]} function.
#'
#' @seealso The diversiting partitioning function \code{\link{partition}}.
#'
#' @import ggplot2
#' @importFrom methods is
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' ## Plot a stacked bar chart of the additive partitioned data
#' partiplot(part.obj, beta.type = "add", plot.type = "bar")
#'
#' ## Plot a line plot of the multiplicative partitioned data
#' partiplot(part.obj, beta.type = "mult", plot.type = "line")
#' }
#' @export
#'
partiplot <- function(part.obj, beta.type = "add", plot.type = "bar"){
  if(!is(part.obj,"partition")){
    stop("'part.obj' must be the set output of the 'partition' function")
  }
  if(!is.character(beta.type)){
    stop("'beta.type' must be set to either 'add' or 'mult'")
  }
  if(!is.character(plot.type)){
    stop("'plot.type' must be set to either 'bar' or 'line'")
  }

  add.beta.rand <- NULL; mult.beta.rand <- NULL; add.beta.rand.ci <- NULL
  mult.beta.rand.ci <- NULL
  levels <- (length(part.obj$Div)-2)/2
  alpha <- as.numeric(part.obj$Div[2])
  add.beta <- as.numeric(part.obj$Div[3:(levels+2)])
  mult.beta <- as.numeric(part.obj$Div[(levels+3):((levels*2)+2)])
  rand.alpha <- mean(as.numeric(part.obj$Rand.Alpha))
  rand.alpha.ci <- 1.96*sd(as.numeric(part.obj$Rand.Alpha))
  for(i in 1:(levels)){
    if(part.obj$Test == "SAMPLE"){
      add.beta.rand[i] <- mean(as.numeric(part.obj$Rand.Beta.Add[[i]]))
      add.beta.rand.ci[i] <- 1.96*sd(as.numeric(part.obj$Rand.Beta.Add[[i]]))
    } else {
      add.beta.rand[i] <- mean(as.numeric(part.obj$Rand.Beta.Add[i,]))
      add.beta.rand.ci[i] <- 1.96*sd(as.numeric(part.obj$Rand.Beta.Add[i,]))

    }
    mult.beta.rand[i] <- mean(part.obj$Rand.Beta.Mult[[i]])
    mult.beta.rand.ci[i] <- 1.96*sd(part.obj$Rand.Beta.Mult[[i]])
  }
  beta.vec <- c(alpha, add.beta, rand.alpha, add.beta.rand, alpha, mult.beta, rand.alpha, mult.beta.rand)
  beta.vec.cu <- c(alpha, add.beta, (rand.alpha + rand.alpha.ci),
                   (add.beta.rand + add.beta.rand.ci), alpha, mult.beta,
                   (rand.alpha + rand.alpha.ci),
                   (mult.beta.rand + mult.beta.rand.ci))
  beta.vec.cl <- c(alpha, add.beta, (rand.alpha - rand.alpha.ci),
                   (add.beta.rand - add.beta.rand.ci), alpha, mult.beta,
                   (rand.alpha - rand.alpha.ci),
                   (mult.beta.rand - mult.beta.rand.ci))
  beta.st <- c(paste(c(rep("Beta",(levels))),1:(levels),sep=" "))
  multidat <- data.frame("DiversityL" = rep(c("Alpha",beta.st,"Alpha",beta.st),times=2),
                         "Diversity"=beta.vec,
                         "Conf.U" = beta.vec.cu,
                         "Conf.L" = beta.vec.cl,
                         "Type" = rep(c("Additive","Multiplicative"),each=((levels*2)+2)),
                         "OvE" = rep(c("Observed","Expected"),each=(levels+1),times=2))
  multidat$OvE <- factor(multidat$OvE,levels(multidat$OvE)[c(2,1)])

  if(plot.type == "bar"){
    if(beta.type == "add"){
      add.bar.plot(multidat)
    } else if(beta.type == "mult"){
      mult.bar.plot(multidat)
    }
  } else if(plot.type == "line"){
    if(beta.type == "add"){
      add.line.plot(multidat)
    } else if(beta.type == "mult"){
      mult.line.plot(multidat)
    }
  }
}

##
add.bar.plot <- function(multidat){
  addbeta <- subset(multidat ,multidat$Type == "Additive")
  addbeta$DiversityL <- factor(addbeta$DiversityL,levels=addbeta$DiversityL[nlevels(addbeta$DiversityL):1])
  addbeta$OvE <- factor(addbeta$OvE, levels = c("Observed", "Expected"))
  ggplot(data = addbeta,aes(x=addbeta$OvE,y=addbeta$Diversity,fill=addbeta$DiversityL))+
    geom_bar(position="fill",stat="identity",color="black")+
    scale_y_continuous(expand=c(0,0))+
    ylab("Total Proportional Diversity")+
    scale_fill_grey(start = 1, end = .2)+
    guides(fill = guide_legend(override.aes = list(size=5)))+
    theme(
      plot.title=element_text(size=rel(1.5)),panel.background=element_rect(
        fill="white",color="black"), axis.text=element_text(size= 18, color="black"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      strip.text.x = element_text(size=16),
      strip.background = element_rect(color="white", fill="white"),
      axis.title=element_text(size=24),
      legend.title=element_text(size=14),
      legend.text=element_text(size=14),
      legend.key=element_blank())
}

##
mult.bar.plot <- function(multidat){
  multbeta <- subset(multidat ,multidat$Type == "Multiplicative")
  multbeta$DiversityL <- factor(multbeta$DiversityL,levels=multbeta$DiversityL[nlevels(multbeta$DiversityL):1])
  multbeta$OvE <- factor(multbeta$OvE, levels = c("Observed", "Expected"))
  ggplot(data = multbeta, aes(x=multbeta$OvE,y=multbeta$Diversity,fill=multbeta$DiversityL))+
    geom_bar(position="fill",stat="identity",color="black")+
    scale_y_continuous(expand=c(0,0))+
    ylab("Total Proportional Diversity")+
    scale_fill_grey(start = 1, end = .2)+
    guides(fill = guide_legend(override.aes = list(size=5)))+
    theme(
      plot.title=element_text(size=rel(1.5)),panel.background=element_rect(
        fill="white",color="black"), axis.text=element_text(size= 18, color="black"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      strip.text.x = element_text(size=16),
      strip.background = element_rect(color="white", fill="white"),
      axis.title=element_text(size=24),
      legend.title=element_text(size=14),
      legend.text=element_text(size=14),
      legend.key=element_blank())
}

##
add.line.plot <- function(multidat){
  addbeta <- subset(multidat ,multidat$Type == "Additive")
  ggplot(data = addbeta,aes(x = addbeta$DiversityL,
                        y = addbeta$Diversity,
                        group = addbeta$OvE)) +
    geom_point(aes(shape = addbeta$OvE), size = 3, fill = "black") +
    stat_summary(fun.y = sum, geom = "line", aes(linetype = addbeta$OvE),
                 size = 1) +
    geom_errorbar(aes(ymin = addbeta$Conf.L,
                      ymax = addbeta$Conf.U,
                      alpha=addbeta$OvE),
                  width = .1) +
    scale_shape_manual(values = c(21, 5)) +
    scale_alpha_manual(values = c(0, 1), guide = FALSE)+
    expand_limits(y = 0) +
    ylab("Diversity")+
    theme(panel.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = 18, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 16),
          strip.background = element_rect(color = "white", fill = "white"),
          axis.title = element_text(size = 24),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key = element_blank())
}

##
mult.line.plot <- function(multidat){
  multbeta <- subset(multidat ,multidat$Type == "Multiplicative")
  ggplot(multbeta, aes(x = multbeta$DiversityL,
                       y = multbeta$Diversity,
                       group = multbeta$OvE)) +
    geom_point(aes(shape = multbeta$OvE), size = 3, fill = "black") +
    stat_summary(fun.y = sum, geom = "line", aes(linetype = multbeta$OvE),
                 size = 1) +
    geom_errorbar(aes(ymin = multbeta$Conf.L,
                      ymax = multbeta$Conf.U,
                      alpha=multbeta$OvE),
                  width = .1) +
    scale_shape_manual(values = c(21, 5)) +
    scale_alpha_manual(values = c(0, 1), guide = FALSE)+
    expand_limits(y = 0) +
    ylab("Diversity")+
    theme(panel.background = element_rect(fill = "white", color = "black"),
          axis.text = element_text(size = 18, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          strip.text.x = element_text(size = 16),
          strip.background = element_rect(color = "white", fill = "white"),
          axis.title = element_text(size = 24),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key = element_blank())
}
