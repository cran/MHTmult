#' Adjusted P-Values for the Double FDR Procedure
#'
#' Given a list/data frame of grouped p-values, retruns adjusted p-values to make decisions
#'
#' @usage
#' DFDR.p.adjust(pval, t, make.decision, alpha)
#' @param pval the structural p-values, the type should be \code{"list"}.
#' @param t the threshold selecting significant families.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}.
#' @return
#' A list of the adjusted p-values, a list of \code{NULL} means the family is not selected to do the test in the second stage.
#' @seealso  \code{\link{DFDR2.p.adjust}}, \code{\link[stats]{p.adjust}}.
#' @author Yalin Zhu
#' @references Mehrotra, D. V., & Heyse, J. F. (2004).
#' Use of the false discovery rate for evaluating clinical safety data.
#' \emph{Statistical methods in medical research}, \strong{13}: 227-238.
#' @examples
#' # data is from Example 4.1 in Mehrotra and Adewale (2012)
#' pval <- list(c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077),
#'              c(0.216,0.843,0.864),
#'              c(1,0.878,0.766,0.598,0.011,0.864),
#'              c(0.889,0.557,0.767,0.009,0.644),
#'              c(1,0.583,0.147,0.789,0.217,1,0.02,0.784,0.579,0.439),
#'              c(0.898,0.619,0.193,0.806,0.611,0.526,0.702,0.196))
#' DFDR.p.adjust(pval = pval,t=0.1)
#' sum(unlist(DFDR.p.adjust(pval = pval,t=0.1))<=0.1)
#' @export
DFDR.p.adjust <- function(pval, t, make.decision=FALSE, alpha=0.05){
  # if (class(pval)!="data.frame" & class(pval)!="list"){
  #   stop("The class of pval has to be 'data.frame' or 'list'.")
  # }
  # if (class(pval)=="data.frame"){
  #   colnames(pval) <- c("family","individual")
  #   pval <- dlply(pval, .(family))
  # }
  adjcondP <- vector("list", length(pval))
  for (i in which(p.adjust(sapply(pval,min),method="BH")<=t) ){
    adjcondP[[i]] <- p.adjust(pval[[i]], method = "BH")
  }
  if (make.decision==FALSE){return(adjcondP)}
  else if(make.decision==TRUE){
    f <- function(x){ifelse(x<=alpha,"reject","not reject")}
    return(list(adjust.p = adjcondP,sig.level=alpha,dec.make=lapply(adjcondP, f)))
  }
}

#' Adjusted P-Values for the Modified Double FDR Procedure
#'
#' Given a list/data frame of grouped p-values, retruns adjusted p-values to make decisions
#'
#' @usage
#' DFDR2.p.adjust(pval, t, make.decision)
#' @param pval the structural p-values, the type should be \code{"list"}.
#' @param t the threshold selecting significant families and testing hypotheses.
#' @param make.decision logical; if \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}.
#' @return
#' A list of the adjusted p-values, a list of \code{NULL} means the family is not selected to do the test in the second stage.
#' @seealso  \code{\link{DFDR.p.adjust}}, \code{\link[stats]{p.adjust}}.
#' @author Yalin Zhu
#' @references Mehrotra, D. V., & Adewale, A. J. (2012).
#' Flagging clinical adverse experiences: reducing false discoveries without materially compromising power for detecting true signals.
#' \emph{Statistics in medicine}, \strong{31}: 1918-1930.
#' @examples
#' # data is from Example 4.1 in Mehrotra and Adewale (2012)
#' pval <- list(c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077),
#'              c(0.216,0.843,0.864),
#'              c(1,0.878,0.766,0.598,0.011,0.864),
#'              c(0.889,0.557,0.767,0.009,0.644),
#'              c(1,0.583,0.147,0.789,0.217,1,0.02,0.784,0.579,0.439),
#'              c(0.898,0.619,0.193,0.806,0.611,0.526,0.702,0.196))
#' DFDR2.p.adjust(pval = pval,t=0.1)
#' sum(unlist(DFDR2.p.adjust(pval = pval,t=0.1))<=0.1)
#' @export
DFDR2.p.adjust <- function(pval, t, make.decision=FALSE){
  # if (class(pval)!="data.frame" & class(pval)!="list"){
  #   stop("The class of pval has to be 'data.frame' or 'list'.")
  # }
  # if (class(pval)=="data.frame"){
  #   colnames(pval) <- c("family","individual")
  #   pval <- dlply(pval, .(family))
  # }
  adjcondP <- vector("list", length(pval))
  adj.pval <- lapply(pval, FUN = p.adjust, method="BH")
  for (i in which(p.adjust(sapply(adj.pval,min),method="BH")<=t) ){
    adjcondP[[i]] <- p.adjust(pval[[i]], method = "BH")
  }
  if (make.decision==FALSE){return(adjcondP)}
  else if(make.decision==TRUE){
    f <- function(x){ifelse(x<=t,"reject","not reject")}
    return(list(adjust.p = adjcondP,sig.level=t,dec.make=lapply(adjcondP, f)))
  }
}
