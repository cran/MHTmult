% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kFWER_adjP.R
\name{gbonf.p.adjust}
\alias{gbonf.p.adjust}
\title{Adjusted P-Values for the Generalized Bonferroni Procedure Controlling k-FWER}
\usage{
gbonf.p.adjust(p, k, alpha, make.decision)
}
\arguments{
\item{p}{numeric vector of p-values (possibly with \code{\link[base]{NA}}s). Any other R is coerced by \code{\link[base]{as.numeric}}. Same as in \code{\link[stats]{p.adjust}}.}

\item{k}{number of allowed type 1 errors in k-FWER controls.}

\item{alpha}{significant level used to compare with adjusted p-values to make decisions, the default value is 0.05.}

\item{make.decision}{logical; if  \code{TRUE}, then the output include the decision rules compared adjusted p-values with significant level \eqn{\alpha}}
}
\value{
A numeric vector of the adjusted p-values (of the same length as \code{p}) if  \code{make.decision = FALSE}, or a list including original p-values, adjusted p-values and decision rules if \code{make.decision = TRUE}.
}
\description{
The function for computing the adjusted p-values based on original p-values and fold \eqn{k}.
}
\examples{
p <- c(0.031,0.023,0.029,0.005,0.031,0.000,0.874,0.399,0.293,0.077)
gbonf.p.adjust(p, k=2)
}
\references{
Lehmann, E. L., & Romano, J. P. (2005).
Generalizations of the familywise error rate.
\emph{The Annals of Statistics}, \strong{33}: 1138-1154.
}
\seealso{
\code{\link{gsidak.p.adjust}}, \code{\link[stats]{p.adjust}}, \code{\link[MHTdiscrete]{Sidak.p.adjust}}.
}
\author{
Yalin Zhu
}
