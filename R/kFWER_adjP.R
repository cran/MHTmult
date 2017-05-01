#' Critical Value for the generalized Sidak Procedure Controlling k-FWER
#'
#' The function for computing the critical value based on number of hypotheses \eqn{m}, fold \eqn{k} and significant level \eqn{\alpha}.
#'
#' @usage
#' gsidak.cv(m, k, alpha)
#' @param m number of hypotheses to be tested.
#' @param k number of allowed type 1 errors in k-FWER controls.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso  \code{\link{gsidak.p.adjust}}, \code{\link[stats]{p.adjust}}, \code{\link[MHTdiscrete]{Sidak.p.adjust}}.
#' @author Yalin Zhu
#' @examples
#' p <- c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077)
#' gsidak.cv(m=length(p), k=2)
#' @export
gsidak.cv <- function(m,k,alpha=0.05) {
  fun <- function(x){
    1-pbinom(q = k-1, size = m, prob = x) - alpha
  }
  uniroot(fun,c(0,1))$root
}

#' Critical Value for the generalized Bonferroni Procedure Controlling k-FWER
#'
#' The function for computing the critical value based on number of hypotheses \eqn{m}, fold \eqn{k} and significant level \eqn{\alpha}.
#'
#' @usage
#' gbonf.cv(m, k, alpha)
#' @param m number of hypotheses to be tested.
#' @param k number of allowed type 1 errors in k-FWER controls.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso  \code{\link{gbonf.p.adjust}}, \code{\link[stats]{p.adjust}}, \code{\link[MHTdiscrete]{Sidak.p.adjust}}.
#' @author Yalin Zhu
#' @examples
#' p <- c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077)
#' gbonf.cv(m=length(p), k=2)
#' @export
gbonf.cv <-  function(m,k,alpha=0.05) {
  k*alpha/m
}

#' Adjusted P-Values for the Generalized Sidak Procedure Controlling k-FWER
#'
#' The function for computing the adjusted p-values based on original p-values and fold \eqn{k}.
#'
#' @usage
#' gsidak.p.adjust(p, k, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param k number of allowed type 1 errors in k-FWER controls.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso  \code{\link{gbonf.p.adjust}}, \code{\link[stats]{p.adjust}}, \code{\link[MHTdiscrete]{Sidak.p.adjust}}.
#' @author Yalin Zhu
#' @references Guo, W., & Romano, J. (2007).
#' A generalized Sidak-Holm procedure and control of generalized error rates under independence.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, \strong{6}(1).
#' @examples
#' p <- c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077)
#' gsidak.p.adjust(p, k=2)
#' @export
gsidak.p.adjust <- function(p, k=1, alpha=0.05, make.decision=FALSE){
  adjP <- 1-pbinom(q = k-1, size = length(p), prob = p)
  if (make.decision==FALSE){
    return(adjP)
  } else{
    return(list(method= paste(k,"-Sidak",sep = ""), significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP, decision=ifelse(adjP<=alpha, "reject","accept"))))
  }
  }

#' Adjusted P-Values for the Generalized Bonferroni Procedure Controlling k-FWER
#'
#' The function for computing the adjusted p-values based on original p-values and fold \eqn{k}.
#'
#' @usage
#' gbonf.p.adjust(p, k, alpha, make.decision)
#' @param p numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.
#' @param k number of allowed type 1 errors in k-FWER controls.
#' @param alpha significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.
#' @param make.decision logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}
#' @return
#' A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
#' @seealso  \code{\link{gsidak.p.adjust}}, \code{\link[stats]{p.adjust}}, \code{\link[MHTdiscrete]{Sidak.p.adjust}}.
#' @author Yalin Zhu
#' @references Lehmann, E. L., & Romano, J. P. (2005).
#' Generalizations of the familywise error rate.
#' \emph{The Annals of Statistics}, \strong{33}: 1138-1154.
#' @examples
#' p <- c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077)
#' gbonf.p.adjust(p, k=2)
#' @export
gbonf.p.adjust <- function(p, k=1, alpha=0.05, make.decision=FALSE){
  adjP <- pmin(length(p)*p/k,1)
  if (make.decision==FALSE){
    return(adjP)
  } else{
    return(list(method= paste(k,"-Bonferroni",sep = ""), significant.level = alpha, Result = data.frame(raw.p = p, adjust.p=adjP, decision=ifelse(adjP<=alpha, "reject","accept"))))
  }
  }
