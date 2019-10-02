\name{BSOIFull}
\alias{BSOIFull}
\title{Bayesian scalar on image (BSOI) regression with fully observed response}
\description{
 This function conducts Bayesian inference for a scalar on image regression with fully observed response.
}
\usage{
 BSOIFull(Y, xi, Z_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0,
 gamma_q0, sigma2_gammaq0, a_delta0, b_delta0, alpha, 
 beta, gamma_cov, sigma2_delta)
}
\arguments{
  \item{Y}{a vector of response with sample size n}
  \item{xi}{a n*K matrix of eigenscores from FPCA}
  \item{Z_cov}{a n*q matrix of covariates}
  \item{n_iter}{number of iterations of the MCMC procedure, default = 10000}
  \item{alpha_0}{mean hyperparameter for alpha, default = 0}
  \item{sigma2_alpha0}{variance hyperparameter for alpha, default = 100}
  \item{beta_k0}{vector of mean hyperparameters for beta, default = rep(0,K)}
  \item{sigma2_betak0}{vector of variance hyperparameters for beta, default = rep(100,K)}
  \item{gamma_q0}{vector of mean hyperparameters for gamma, default = rep(0.01,q)}
  \item{sigma2_gammaq0}{vector of variance hyperparameters for gamma, default = rep(100,q)}
  \item{a_delta0}{hyperparameters for the error variance, default = 9}
  \item{b_delta0}{hyperparameters for the error variance, default = 3}
  \item{alpha, beta, gamma_cov, sigma2_delta}{initial values of the parameters, with default values set as, 0.1, rep(0.1,K), rep(0.1,q), and 1, respectively.}
}
\value{
  \code{BSOIFull()} returns a list that contains posterior samples for the parameters, i.e., alpha, beta, gamma_cov, and sigma2_delta.
}

\references{
  Feng X.N., Li T.F., Song X.Y. and Zhu H.T. (2017) Bayesian Scalar on Image Regression with Non-ignorable Non-response.
}

\examples{
#This example applies default values for the prior inputs and initial values
#BSOIFull(Y, xi, Z_cov)

##please also see the democode for the simulation study
}


