#' Selecting Threshold for cFDR Controlling Procedures
#'
#' Given the structural p-values, choose a selecting method for controlling generalized familywise error rate or false discovery rate across families, and a combining mehtod, returns a vector of thresholds for the first stage of cFDR controlling procedures.
#'
#' @usage
#' select.thres(pval, select.method, comb.method, alpha, k)
#' @param pval the structural p-values, the type should be \code{"list"}.
#' @param select.method global p-value selecting methods. For generalized FWER controlling, k-Bonferroni or k-Sidak procedures can be used; for FDR controlling, BH procedure can be used.
#' @param comb.method p-value combining methods including \code{"Fisher"}, \code{"Stouffer"}, and \code{"minP"} combining methods.
#' @param alpha significant level for selecting significant families in the first stage. The default value is 0.05.
#' @param k number of allowed type 1 errors in k-FWER controls.
#' @return A list of the adjusted conditional p-values, a list of \code{NULL} means the family is not selected to do the test in the second stage.
#' @author Yalin Zhu
#' @export
select.thres <- function(pval, select.method=c("bonferroni","sidak","BH"),
                         comb.method=c("Fisher","Stouffer","minP"), alpha=0.05, k=1){
  select.method <- match.arg(select.method)
  comb.method <- match.arg(comb.method)
  # if (class(pval)!="data.frame" & class(pval)!="list"){
  #   stop("The class of pval has to be 'data.frame' or 'list'.")
  # }
  # if (class(pval)=="data.frame"){
  #   colnames(pval) <- c("family","individual")
  #   pval <- dlply(pval, .(family))
  # }
   if (comb.method=="minP"){
    global.p <- sapply(pval, length)*sapply(pval,min)
    if (select.method=="sidak"){
      t <- 1-(1-gsidak.cv(m = length(global.p),k = k,alpha = alpha))^(1/sapply(pval,length))
    }
    else if (select.method=="bonferroni"){
  #    t <- k*alpha/(length(global.p)*sapply(pval,length))
      t <- 1-(1-k*alpha/length(global.p))^(1/sapply(pval,length))
    }
    else if (select.method=="BH"){
      t <- c()
      for (i in 1:length(global.p)){
        J <- max(sum(p.adjust(global.p[-i],method="BH") <= alpha),1)
  #      t[i] <- (J+1)/length(global.p)*alpha/length(pval[[i]])
        t[i] <- 1- (1-(J+1)/length(global.p)*alpha)^(1/length(pval[[i]]))
      }
    }
  }
  # Fisher's combining method
  else if (comb.method=="Fisher"){
  f <- function(x){sum(log(x))}

  global.p <- 1-pchisq(-2*sapply(pval,FUN = f), 2*sapply(pval,length))

  if (select.method=="sidak"){
    t <- qchisq(1-gsidak.cv(m = length(global.p),k = k,alpha = alpha), 2*sapply(pval,length))
  }

  else if (select.method=="bonferroni"){
    t <- qchisq(1-k*alpha/length(global.p), 2*sapply(pval,length))
  }
  else if (select.method=="BH"){
  t <- c()
    for (i in 1:length(global.p)){
      J <- sum(p.adjust(global.p[-i],method="BH") <= alpha)
      t[i] <- qchisq(1-(J+1)/length(global.p)*alpha, 2*length(pval[[i]]))
    }
  }
}
  return(t)

}


#' Adjusted Conditional P-values for Two-stage cFDR Controlling Procedures
#'
#' Given a list/data frame of grouped p-values, selecting thresholds and p-value combining method, retruns adjusted conditional p-values to make decisions
#'
#' @usage
#' cFDR.cp.adjust(pval, t, comb.method = c("Fisher", "Stouffer", "minP"),
#' make.decision, sig.level)
#' @param pval the structural p-values, the type should be \code{"list"}.
#' @param t the thresholds determining whether the families are selected or not, also affects conditional p-value within families.
#' @param comb.method p-value combining methods including \code{"Fisher"}, \code{"Stouffer"}, and \code{"minP"} combining methods.
#' @param make.decision logical; if \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}.
#' @param sig.level significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @import stats
#' @return A list of the adjusted conditional p-values, a list of \code{NULL} means the family is not selected to do the test in the second stage.
#' @author Yalin Zhu
#' @references Heller, R., Chatterjee, N., Krieger, A., & Shi, J. (2016).
#' Post-selection Inference Following Aggregate Level Hypothesis Testing in Large Scale Genomic Data.
#' \emph{bioRxiv}, 058404.
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
#' t=select.thres(pval,select.method = "BH", comb.method = "minP", alpha = 0.1)
#' cFDR.cp.adjust(pval, t=t, comb.method="minP")
#'
#' t1=select.thres(pval, select.method = "bonferroni", comb.method = "minP", alpha = 0.1, k=3)
#' cFDR.cp.adjust(pval, t=t1, comb.method="minP")
#'
#' t2=select.thres(pval, select.method = "sidak", comb.method = "minP", alpha = 0.1, k=3)
#' cFDR.cp.adjust(pval, t=t2, comb.method="minP")
#' @export
cFDR.cp.adjust <- function(pval, t, comb.method=c("Fisher", "Stouffer", "minP"),  make.decision=FALSE, sig.level=0.05){
  comb.method <- match.arg(comb.method)
  if (length(t)<length(pval)){t <- rep(t,len=length(pval))}
  condPval <- vector("list", length(pval))
  adjcondP <- vector("list", length(pval))

  if (comb.method=="Fisher"){
    for (i in which(lapply(pval,prod)<=exp(-t/2)) ){
      condPval[[i]] <- integer(length(pval[[i]]))
      for (j in 1:length(pval[[i]])){
        if (length(pval[[i]])==1){condPval[[i]][j] <- pval[[i]][j]/exp(-t[i]/2)}
        else {
          condPval[[i]][j] <- pval[[i]][j]/min(c(exp(-t[i]/2)/prod(pval[[i]][-j]),1))
        }
      }
      adjcondP[[i]] <- p.adjust(condPval[[i]], method = "BH")
    }
  }

  if (comb.method=="Stouffer"){
    stfun <- function(x){sum(qnorm(1-x))/sqrt(length(x))}
    for (i in which(sapply(pval,stfun)>=t) ){
      condPval[[i]] <- integer(length(pval[[i]]))
      for (j in 1:length(pval[[i]])){
        if (length(pval[[i]])==1){condPval[[i]][j] <- pval[[i]][j]/(1-pnorm(t[i]))}
        else {
        condPval[[i]][j] <- pval[[i]][j]/(1-pnorm(sqrt(length(pval[[i]]))*t[i]-sum(qnorm(1-pval[[i]][-j]))))
        }
      }
      adjcondP[[i]] <- p.adjust(condPval[[i]], method = "BH")
    }
  }

  else if (comb.method=="minP"){
    for (i in which(sapply(pval,min)<=t) ){
      condPval[[i]] <- integer(length(pval[[i]]))
      for (j in 1:length(pval[[i]])){
        if (length(pval[[i]])==1){condPval[[i]][j] <- pval[[i]][j]/t[i]}
        else {
          condPval[[i]][j] <- pval[[i]][j]/ifelse(min(pval[[i]][-j]) > t[i], t[i] ,1)
        }
      }
      adjcondP[[i]] <- p.adjust(condPval[[i]], method = "BH")
    }
  }
  if (make.decision==FALSE){return(adjcondP)}
  else if(make.decision==TRUE){
    f <- function(x){ifelse(x<=sig.level,"reject","not reject")}
    return(list(adjust.p = adjcondP,adjust.sig.level=sig.level,dec.make=lapply(adjcondP, f)))
  }
}
