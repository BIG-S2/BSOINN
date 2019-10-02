set.seed(as.numeric(Sys.time()))
library(BSOINN)  ##load the library

####file names
Y_file = "Y.txt"
Cov_file = "COV.txt"  
R_file = "R.txt"
Xi_file = "XI.txt"

####define constants
N_REPS = 1
N_subj=802     ##sample size
n_iter=20000    ##total iterations
n_burnin=10000  ##burn-in phase
K_x = 5        ##number of eigenimages retained   
q_cov = 7      ##dimension of covariates in SoI
q_cov_NN = 6   ##dimension of covariates in missing mechanism after removing the instrument 
ind_cov_NN = c(T,T,F,T,T,T,T)  ##we define educational level (the third) as an instrument


####define objects to record the estimation results
Ealpha_IN = array(0, dim = c(N_REPS,1))
Egamma_IN = array(0,dim = c(N_REPS,q_cov))
Ebeta_IN = array(0, dim = c(N_REPS,K_x))
Edelta_IN = array(0, dim = c(N_REPS,1))

SDalpha_IN = array(0, dim = c(N_REPS,1))
SDgamma_IN = array(0,dim = c(N_REPS,q_cov))
SDbeta_IN = array(0, dim = c(N_REPS,K_x))
SDdelta_IN = array(0, dim = c(N_REPS,1))

Ealpha_NN = array(0, dim = c(N_REPS,1))
Egamma_NN = array(0,dim = c(N_REPS,q_cov))
Ebeta_NN = array(0, dim = c(N_REPS,K_x))
Edelta_NN = array(0, dim = c(N_REPS,1))
Ealpha_r_NN = array(0, dim = c(N_REPS,1))
Egamma_r_NN = array(0,dim = c(N_REPS,q_cov_NN))
Ebeta_r_NN =array(0, dim = c(N_REPS,K_x))
Ephi_NN = array(0, dim = c(N_REPS,1))

SDalpha_NN = array(0, dim = c(N_REPS,1))
SDgamma_NN = array(0,dim = c(N_REPS,q_cov))
SDbeta_NN = array(0, dim = c(N_REPS,K_x))
SDdelta_NN = array(0, dim = c(N_REPS,1))
SDalpha_r_NN = array(0, dim = c(N_REPS,1))
SDgamma_r_NN = array(0,dim = c(N_REPS,q_cov_NN))
SDbeta_r_NN =array(0, dim = c(N_REPS,K_x))
SDphi_NN = array(0, dim = c(N_REPS,1))

varphi_acceptrate_NN = array(0, dim = c(N_REPS,1))
ym_acceptrate_NN = array(0,dim=c(N_REPS, 1))

####analyze the real data
for(CIR in 1:N_REPS){
	####read the values
	Y = matrix(scan(Y_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1, byrow = T)   ##outcome subject to missingness
	X_cov = matrix(scan(Cov_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=q_cov, byrow = T)  ##covariates in SoI regression
	R = matrix(scan(R_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1,byrow = T)  ##missing indicator
	xi = matrix(scan(Xi_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=K_x, byrow=T)  ##eigenscores from FPCA, and the FPCA could be easily conducted through matlab or the fast.svd function in R
	cat("\n")
	print(paste("Replication:", CIR))
	
	X_cov_NN = X_cov[,ind_cov_NN]
	N_mis = sum(R==1)
	
	####SoI_IN procedure
	print("BSOI_IN:")
	##hyperparameters
	alpha_0 = 0
	sigma2_alpha0 = 100
	beta_k0 = rep(0,K_x)
	sigma2_betak0 = rep(100,K_x)
	gamma_q0 = rep(0.01,q_cov)
	sigma2_gammaq0 = rep(100,q_cov)
	a_delta0 = 9
	b_delta0 = 3
    
	##temp values for missing components
	Y[R==1] = rnorm(N_mis,mean = 0, sd = 1)
	
	##initial values for the patameters, set 1
	alpha_set1 = 0.1
	beta_set1 = rep(0.1,K_x)
	gamma_cov_set1 = rep(0.1,q_cov)
	sigma2_delta_set1 = 1
	
	
	##initial values for the parameters, set 2
	alpha_set2 = 1
	beta_set2 = rep(1,K_x)
	gamma_cov_set2 = rep(1,q_cov)
	sigma2_delta_set2 = 1
	
	##initial values for the parameters, set 3
	alpha_set3 = -1
	beta_set3 = rep(-1,K_x)
	gamma_cov_set3 = rep(-1,q_cov)
	sigma2_delta_set3 = 1
	
	##conduct the Bayesian analysis using the function SoIIN()
	##Only the data inputs are compulsory, and we may implement the analysis as
    ##model_IN <- BSOIIN(Y, R, xi, X_cov, n_iter)
	##The convergence of the model is evaluated via impose different initial values on the parameters.    
	model_IN1 <- BSOIIN(Y, R, xi, X_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0,   gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set1, beta_set1,  gamma_cov_set1, sigma2_delta_set1)
	model_IN2 <- BSOIIN(Y, R, xi, X_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0,   gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set2, beta_set2,  gamma_cov_set2, sigma2_delta_set2)
	model_IN3 <- BSOIIN(Y, R, xi, X_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0,   gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set3, beta_set3,  gamma_cov_set3, sigma2_delta_set3)
	
	#beta and gamma are major parameters, and the convergence is checked through the traceplots of multiple chains
	par(mfrow=c(3,2))
	for(iii in 1:K_x){
		plot(model_IN1$beta[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("beta",iii), main = "BSOIIN" )
		lines(model_IN2$beta[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
		lines(model_IN3$beta[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	par(mfrow=c(4,2))
	for(iii in 1:q_cov){
	  plot(model_IN1$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("gamma",iii), main = "BSOIIN" )
	  lines(model_IN2$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
	  lines(model_IN3$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	##record the results
	Ealpha_IN[CIR,] = mean(c(model_IN1$alpha[(n_burnin+1):n_iter],model_IN2$alpha[(n_burnin+1):n_iter],model_IN3$alpha[(n_burnin+1):n_iter]))
	Ebeta_IN[CIR,] = apply(rbind(model_IN1$beta[(n_burnin+1):n_iter,],model_IN2$beta[(n_burnin+1):n_iter,], model_IN3$beta[(n_burnin+1):n_iter,]),2,mean)
	Egamma_IN[CIR,] = apply(rbind(model_IN1$gamma_cov[(n_burnin+1):n_iter,],model_IN2$gamma_cov[(n_burnin+1):n_iter,],model_IN3$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_IN[CIR,] = mean(c(model_IN1$sigma2_delta[(n_burnin+1):n_iter],model_IN2$sigma2_delta[(n_burnin+1):n_iter],model_IN3$sigma2_delta[(n_burnin+1):n_iter]))
	
	SDalpha_IN[CIR,] = sd(c(model_IN1$alpha[(n_burnin+1):n_iter],model_IN2$alpha[(n_burnin+1):n_iter],model_IN3$alpha[(n_burnin+1):n_iter]))
	SDbeta_IN[CIR,] = apply(rbind(model_IN1$beta[(n_burnin+1):n_iter,],model_IN2$beta[(n_burnin+1):n_iter,], model_IN3$beta[(n_burnin+1):n_iter,]),2,sd)
	SDgamma_IN[CIR,] = apply(rbind(model_IN1$gamma_cov[(n_burnin+1):n_iter,],model_IN2$gamma_cov[(n_burnin+1):n_iter,],model_IN3$gamma_cov[(n_burnin+1):n_iter,]),2,sd)
	SDdelta_IN[CIR,] = sd(c(model_IN1$sigma2_delta[(n_burnin+1):n_iter],model_IN2$sigma2_delta[(n_burnin+1):n_iter],model_IN3$sigma2_delta[(n_burnin+1):n_iter]))
	
	
	##print results
	cat("\n")
	print("#########Estimation results of BSOIIN:")
	print("##beta:")
	print("mean:")
	print(Ebeta_IN[CIR,])
	print("SD:")
	print(SDbeta_IN[CIR,])
	print("##gamma:")
	print("mean:")
	print(Egamma_IN[CIR,])
	print("SD:")
	print(SDgamma_IN[CIR,])

	####SoI_NN procedure
	cat("\n")
	print("BSOI_NN:")
	##variance parameters that control acceptance rate in MH algorithm
	var_ym = 20
	var_varphi = 1.5
	
	##hyperparameters
	alpha_0 = 0
	sigma2_alpha0 = 100
	beta_k0 = rep(0,K_x)
	sigma2_betak0 = rep(100,K_x)
	gamma_q0 = rep(0.01,q_cov)
	sigma2_gammaq0 = rep(100,q_cov)
	a_delta0 = 9
	b_delta0 = 3
	varphi_0 = c(0, rep(0,q_cov_NN), rep(0,K_x), 0)
	sigma_varphi0 = diag(2+q_cov_NN+K_x)
	
	##temp values for missing components
	Y[R==1] = rnorm(N_mis,mean = 0, sd = 1)
	
	##initial values for the parameters
	alpha_set1 = 0.1
	beta_set1 = rep(0.1,K_x)
	gamma_cov_set1 = rep(0.1,q_cov)
	sigma2_delta_set1 = 1
	alpha_r_set1 = 0.1
	beta_r_set1 = rep(0.1,K_x)
	gamma_cov_r_set1 = rep(0.1,q_cov_NN)
	phi_set1 = 0.1
	
	##initial values for the parameters
	alpha_set2 = 1
	beta_set2 = rep(1,K_x)
	gamma_cov_set2 = rep(1,q_cov)
	sigma2_delta_set2 = 1
	alpha_r_set2 = 1
	beta_r_set2 = rep(1,K_x)
	gamma_cov_r_set2 = rep(1,q_cov_NN)
	phi_set2 = 1
	
	##initial values for the parameters
	alpha_set3 = -1
	beta_set3 = rep(-1,K_x)
	gamma_cov_set3 = rep(-1,q_cov)
	sigma2_delta_set3 = 1
	alpha_r_set3 = -1
	beta_r_set3 = rep(-1,K_x)
	gamma_cov_r_set3 = rep(-1,q_cov_NN)
	phi_set3 = -1

	
	##conduct the Bayesian analysis using the function SoINN()
	##Only the data inputs are compulsory, and thus, we may implement the analysis as
    ##model_NN <- BSOINN(Y, R, xi, X_cov,  X_cov_NN, n_iter)
	##The convergence of the model is evaluated via impose different initial values on the parameters.    	
	model_NN1 <- BSOINN(Y, R, xi, X_cov, X_cov_NN, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0, gamma_q0, sigma2_gammaq0, a_delta0, b_delta0, varphi_0, sigma_varphi0,  var_ym,var_varphi, alpha_set1, beta_set1, gamma_cov_set1, sigma2_delta_set1,  alpha_r_set1, beta_r_set1, gamma_cov_r_set1, phi_set1)
	model_NN2 <- BSOINN(Y, R, xi, X_cov, X_cov_NN, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0, gamma_q0, sigma2_gammaq0, a_delta0, b_delta0, varphi_0, sigma_varphi0,  var_ym,var_varphi, alpha_set2, beta_set2, gamma_cov_set2, sigma2_delta_set2,  alpha_r_set2, beta_r_set2, gamma_cov_r_set2, phi_set2)
	model_NN3 <- BSOINN(Y, R, xi, X_cov, X_cov_NN, n_iter, alpha_0, sigma2_alpha0, beta_k0, sigma2_betak0, gamma_q0, sigma2_gammaq0, a_delta0, b_delta0, varphi_0, sigma_varphi0,  var_ym,var_varphi, alpha_set3, beta_set3, gamma_cov_set3, sigma2_delta_set3,  alpha_r_set3, beta_r_set3, gamma_cov_r_set3, phi_set3)
	
	#beta and gamma are major parameters, and the convergence is checked through the traceplots of multiple chains
	par(mfrow=c(3,2))
	for(iii in 1:K_x){
		plot(model_NN1$beta[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("beta",iii), main = "BSOINN" )
		lines(model_NN2$beta[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
		lines(model_NN3$beta[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	par(mfrow=c(4,2))
	for(iii in 1:q_cov){
	  plot(model_NN1$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("gamma",iii), main = "BSOINN" )
	  lines(model_NN2$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
	  lines(model_NN3$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	par(mfrow = c(1,1))
	plot(model_NN1$phi[(n_burnin+1):n_iter],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = "phi", main = "BSOINN", ylim = c(-3,0) )
	lines(model_NN2$phi[(n_burnin+1):n_iter],type = "l",col = "blue")
	lines(model_NN3$phi[(n_burnin+1):n_iter],type = "l",col = "green")
	
	par(mfrow=c(3,2))
	for(iii in 1:K_x){
		plot(model_NN1$beta_r[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("beta_r",iii), main = "BSOINN" )
		lines(model_NN2$beta_r[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
		lines(model_NN3$beta_r[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	par(mfrow=c(4,2))
	for(iii in 1:q_cov_NN){
	  plot(model_NN1$gamma_cov_r[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("gamma_r",iii), main = "BSOINN" )
	  lines(model_NN2$gamma_cov_r[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
	  lines(model_NN3$gamma_cov_r[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	##record results
	Ealpha_NN[CIR,] = mean(c(model_NN1$alpha[(n_burnin+1):n_iter],model_NN2$alpha[(n_burnin+1):n_iter],model_NN3$alpha[(n_burnin+1):n_iter]))
	Ebeta_NN[CIR,] = apply(rbind(model_NN1$beta[(n_burnin+1):n_iter,],model_NN2$beta[(n_burnin+1):n_iter,],model_NN3$beta[(n_burnin+1):n_iter,]),2,mean)
	Egamma_NN[CIR,] = apply(rbind(model_NN1$gamma_cov[(n_burnin+1):n_iter,],model_NN2$gamma_cov[(n_burnin+1):n_iter,],model_NN3$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_NN[CIR,] = mean(c(model_NN1$sigma2_delta[(n_burnin+1):n_iter],model_NN2$sigma2_delta[(n_burnin+1):n_iter],model_NN3$sigma2_delta[(n_burnin+1):n_iter]))
	
	Ealpha_r_NN[CIR,] = mean(c(model_NN1$alpha_r[(n_burnin+1):n_iter],model_NN2$alpha_r[(n_burnin+1):n_iter],model_NN3$alpha_r[(n_burnin+1):n_iter]))
	Ebeta_r_NN[CIR,] = apply(rbind(model_NN1$beta_r[(n_burnin+1):n_iter,],model_NN2$beta_r[(n_burnin+1):n_iter,],model_NN3$beta_r[(n_burnin+1):n_iter,]),2,mean)
	Egamma_r_NN[CIR,] = apply(rbind(model_NN1$gamma_cov_r[(n_burnin+1):n_iter,],model_NN2$gamma_cov_r[(n_burnin+1):n_iter,],model_NN3$gamma_cov_r[(n_burnin+1):n_iter,]),2,mean)
	Ephi_NN[CIR,] = mean(c(model_NN1$phi[(n_burnin+1):n_iter],model_NN2$phi[(n_burnin+1):n_iter],model_NN3$phi[(n_burnin+1):n_iter]))
	
	SDalpha_NN[CIR,] = sd(c(model_NN1$alpha[(n_burnin+1):n_iter],model_NN2$alpha[(n_burnin+1):n_iter],model_NN3$alpha[(n_burnin+1):n_iter]))
	SDbeta_NN[CIR,] = apply(rbind(model_NN1$beta[(n_burnin+1):n_iter,],model_NN2$beta[(n_burnin+1):n_iter,],model_NN3$beta[(n_burnin+1):n_iter,]),2,sd)
	SDgamma_NN[CIR,] = apply(rbind(model_NN1$gamma_cov[(n_burnin+1):n_iter,],model_NN2$gamma_cov[(n_burnin+1):n_iter,],model_NN3$gamma_cov[(n_burnin+1):n_iter,]),2,sd)
	SDdelta_NN[CIR,] = sd(c(model_NN1$sigma2_delta[(n_burnin+1):n_iter],model_NN2$sigma2_delta[(n_burnin+1):n_iter],model_NN3$sigma2_delta[(n_burnin+1):n_iter]))
	
	SDalpha_r_NN[CIR,] = sd(c(model_NN1$alpha_r[(n_burnin+1):n_iter],model_NN2$alpha_r[(n_burnin+1):n_iter],model_NN3$alpha_r[(n_burnin+1):n_iter]))
	SDbeta_r_NN[CIR,] = apply(rbind(model_NN1$beta_r[(n_burnin+1):n_iter,],model_NN2$beta_r[(n_burnin+1):n_iter,],model_NN3$beta_r[(n_burnin+1):n_iter,]),2,sd)
	SDgamma_r_NN[CIR,] = apply(rbind(model_NN1$gamma_cov_r[(n_burnin+1):n_iter,],model_NN2$gamma_cov_r[(n_burnin+1):n_iter,],model_NN3$gamma_cov_r[(n_burnin+1):n_iter,]),2,sd)
	SDphi_NN[CIR,] = sd(c(model_NN1$phi[(n_burnin+1):n_iter],model_NN2$phi[(n_burnin+1):n_iter],model_NN3$phi[(n_burnin+1):n_iter]))
	
	#acceptance rate in the MH algorithm
	varphi_acceptrate_NN[CIR,] = (model_NN1$acce_varphi_n + model_NN2$acce_varphi_n + model_NN3$acce_varphi_n)/(3*n_iter)
	ym_acceptrate_NN[CIR,] =  (sum(model_NN1$acce_ym)+sum(model_NN2$acce_ym)+sum(model_NN3$acce_ym))/3/n_iter/N_mis
	
	##print results
	cat("\n")
	print("#########Estimation results of BSOINN:")
	print("####BSOI model:")
	print("##beta:")
	print("mean:")
	print(Ebeta_NN[CIR,])
	print("SD:")
	print(SDbeta_NN[CIR,])
	print("##gamma:")
	print("mean:")
	print(Egamma_NN[CIR,])
	print("SD:")
	print(SDgamma_NN[CIR,])
	print("####Tilting model:")
	print("##beta_r:")
	print("mean:")
	print(Ebeta_r_NN[CIR,])
	print("SD:")
	print(SDbeta_r_NN[CIR,])
	print("##gamma_r:")
	print("mean:")
	print(Egamma_r_NN[CIR,])
	print("SD:")
	print(SDgamma_r_NN[CIR,])
	print("##phi:")
	print("mean:")
	print(Ephi_NN[CIR,])
	print("SD:")
	print(SDphi_NN[CIR,])
} 