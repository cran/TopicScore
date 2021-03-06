% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/W_from_AD.R
\name{W_from_AD}
\alias{W_from_AD}
\title{Estimation of W from A and X}
\usage{
W_from_AD(A, X)
}
\arguments{
\item{A}{The p-by-K word-topic matrix.}

\item{X}{The p-by-n word-document matrix.}
}
\value{
The estimated K-by-n topic-document matrix W_hat.
}
\description{
This function estimates the topic-document matrix W from the word-topic matrix A and the word-document X through quadratic programming.
}
\examples{
data("AP")
K <- 3
tscore_obj <- topic_score(K, AP)
W_hat <- W_from_AD(tscore_obj$A_hat, AP)
}
\author{
Minzhe Wang
}
