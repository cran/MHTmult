#' Adjusted P-values for the Group BH Procedure
#'
#' Given a list/data frame of grouped p-values, selecting thresholds and p-value combining method, retruns adjusted conditional p-values to make decisions
#'
#' @usage
#' GBH.p.adjust(pval, t, make.decision)
#' @param pval the structural p-values, the type should be \code{"list"}.
#' @param t the thresholds determining whether the families are selected or not, also affects conditional p-value within families.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return A list of the adjusted conditional p-values, a list of \code{NULL} means the family is not selected to do the test in the second stage.
#' @author Yalin Zhu
#' @references Hu, J. X., Zhao, H., & Zhou, H. H. (2010).
#' False discovery rate control with groups.
#' \emph{Journal of the American Statistical Association}, \strong{105}: 1215-1227.
#' @examples
#' # data is from Example 4.1 in Mehrotra and Adewale (2012)
#' pval <- list(c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077),
#'              c(0.216,0.843,0.864),
#'              c(1,0.878,0.766,0.598,0.011,0.864),
#'              c(0.889,0.557,0.767,0.009,0.644),
#'              c(1,0.583,0.147,0.789,0.217,1,0.02,0.784,0.579,0.439),
#'              c(0.898,0.619,0.193,0.806,0.611,0.526,0.702,0.196))
#' sum(p.adjust(unlist(pval), method = "BH")<=0.1)
#' DFDR.p.adjust(pval = pval,t=0.1)
#' DFDR2.p.adjust(pval = pval,t=0.1)
#' sum(unlist(DFDR.p.adjust(pval = pval,t=0.1))<=0.1)
#' sum(unlist(DFDR2.p.adjust(pval = pval,t=0.1))<=0.1)
#'
#' GBH.p.adjust(pval = pval,t=0.1)
#' sum(unlist(GBH.p.adjust(pval = pval,t=0.1))<=0.1)
#'
#' t=select.thres(pval,select.method = "BH", comb.method = "minP", alpha = 0.1)
#' cFDR.cp.adjust(pval, t=t, comb.method="minP")
#'
#' t1=select.thres(pval, select.method = "bonferroni", comb.method = "minP", alpha = 0.1, k=3)
#' cFDR.cp.adjust(pval, t=t1, comb.method="minP")
#'
#' t2=select.thres(pval, select.method = "sidak", comb.method = "minP", alpha = 0.1, k=3)
#' cFDR.cp.adjust(pval, t=t2, comb.method="minP")
#' @export
GBH.p.adjust <- function(pval, t, make.decision=FALSE){
  # if (class(pval)!="data.frame" & class(pval)!="list"){
  #   stop("The class of pval has to be 'data.frame' or 'list'.")
  # }
  # if (class(pval)=="data.frame"){
  #   colnames(pval) <- c("family","individual")
  #   pval <- dlply(pval, .(family))
  # }
  adj.pval <- lapply(pval, FUN = p.adjust, method="BH")
    f1 <- function(x){sum(x<=t/(1+t))}
  pi0.hat <- 1-sapply(adj.pval, f1)/sapply(pval,length)
    f2 <- function(x,y){pmin(y/(1-y)*x,1)}
  p.weight <- mapply(f2, x=pval, y=pi0.hat)
  adjP <- p.adjust(unlist(p.weight), method = "BH")
    if (make.decision==FALSE){return(adjP)}
    else if(make.decision==TRUE){
    avgpi1.hat <- sum(sapply(adj.pval, f1))/sum(sapply(pval,length))
    sig.level=t/(1+t)/avgpi1.hat
    f <- function(x){ifelse(x<=sig.level,"reject","not reject")}
    return(list(adjust.p = adjP,adjust.sig.level=sig.level,dec.make=f(adjP)))
  }
}


