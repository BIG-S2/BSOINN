\name{BSOINN}
\alias{BSOINN}
\title{Bayesian scalar on image (BSOI) regression with  non-ignorable non-response}
\description{
 This function conducts Bayesian inference for a scalar on image regression with non-ignorable non-response.
}
\usage{
BSOINN(Y, R, xi, Z_cov, Z_cov_NN, n_iter, alpha_0, sigma2_alpha0, beta_k0, 
sigma2_betak0, gamma_q0, sigma2_gammaq0, a_delta0, b_delta0, 
varphi_0, sigma_varphi0, var_ym, var_varphi, 
alpha, beta, gamma_cov, sigma2_delta, alpha_r, beta_r, gamma_cov_r, phi)
}
\arguments{
  \item{Y}{a vector of response with the sample size n}
  \item{R}{an integer vector of missingness indicators (1 indicates missing and 0 indicates observable) for response with the sample size n}
  \item{xi}{a n*K matrix of eigenscores from FPCA}
  \item{Z_cov}{a n*q matrix of covariates}
  \item{Z_cov_NN}{a n*q_mis matrix of covariates (without the instrumental variable) in the exponential tilting model}
  \item{n_iter}{number of iterations of the MCMC procedure, default = 10000}
  \item{alpha_0}{mean hyperparameter for alpha, default = 0}
  \item{sigma2_alpha0}{variance hyperparameter for alpha, default = 100}
  \item{beta_k0}{vector of mean hyperparameters for beta, default = rep(0,K)}
  \item{sigma2_betak0}{vector of variance hyperparameters for beta, default = rep(100,K)}
  \item{gamma_q0}{vector of mean hyperparameters for gamma, default = rep(0.01,q)}
  \item{sigma2_gammaq0}{vector of variance hyperparameters for gamma, default = rep(100,q)}
  \item{a_delta0}{hyperparameters for the error variance, default = 9}
  \item{b_delta0}{hyperparameters for the error variance, default = 3}
  \item{varphi_0, sigma_varphi0}{hyperparameters for the parameters in the missing mechanism, with default values set as, c(0, rep(0,q_mis), rep(0,K), 0), and diag(2+q_mis+K), respectively}
  \item{var_ym, var_varphi}{the variance parameters that control the acceptance rate in the MH algorithm, with default values set as, 20 and 1.5, repectively} 
   \item{alpha, beta, gamma_cov, sigma2_delta}{initial values of the parameters, with default values set as, 0.1, rep(0.1,K), rep(0.1,q), and 1, respectively.}
  \item{alpha_r, beta_r, gamma_cov_r, phi}{initial values of the parameters in the missing mechanism, with default values set as, 0.1, rep(0.1,K), rep(0.1,q_mis), and 0.1, respectively.}
}
\value{
  \code{BSOINN()} returns a list that contains posterior samples for the parameters, i.e., alpha, beta, sigma2_delta, gamma_cov, alpha_r, beta_r, gamma_cov_r, and phi. The list also contains the number of iterations that the candidate value from the proposal distribution is accepted in the MH algorithm for non-response and the parameters in the missing mechanism.
}

\references{
   Feng X.N., Li T.F., Song X.Y. and Zhu H.T. (2017) Bayesian Scalar on Image Regression with Non-ignorable Non-response.
}

\examples{
#This example applies default values for the prior inputs and initial values
#BSOINN(Y, R, xi, Z_cov, Z_cov_NN)

##please also see the democode for the simulation study
}


