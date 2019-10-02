set.seed(as.numeric(Sys.time()))
library(BSOINN)  ##load the library
library(corpcor) ##contains the fast.svd function

####file names
Outcome_file = "OUTCOME.txt" 
Y_file = "Y.txt"
Cov_file = "COV.txt"
R_file = "R.txt"
Xi_file = "XI.txt"
Eigen_ind = "IND_REV.txt"  ##In the simulation study, we need to adjust the sign of the eigenimages 

####define constants
dim_image = 300
N_REPS = 1
N_subj=1000     ##sample size
n_iter=10000    ##total iterations
n_burnin=4000  ##burn-in phase
K_x = 3        ##number of eigenimages retained   
q_cov = 3      ##dimension of covariates in SoI
q_cov_NN = 2   ##dimension of covariates in missing mechanism after removing the instrument 
ind_cov_NN = c(T,T,F)  ##we define the last covariate as the instrument

######True values in data generation process
sigma2_delta_true = 1
alpha_true = 0
beta_true = c(0.5,1,-1)
gamma_true = c(1.5,-1,0.5)

alpha_r_true = 0.5
gamma_r_true = c(-0.7,-0.7)
beta_r_true = c(-1,0.5,0.5)
phi_true = -1.2

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

Ealpha_Full = array(0, dim = c(N_REPS,1))
Egamma_Full = array(0,dim = c(N_REPS,q_cov))
Ebeta_Full = array(0, dim = c(N_REPS,K_x))
Edelta_Full = array(0, dim = c(N_REPS,1))

SDalpha_Full = array(0, dim = c(N_REPS,1))
SDgamma_Full = array(0,dim = c(N_REPS,q_cov))
SDbeta_Full = array(0, dim = c(N_REPS,K_x))
SDdelta_Full = array(0, dim = c(N_REPS,1))

####Data generation
for(CIR in 1:N_REPS){
	####true eigenvectors
	eigen_vector = array(0,c(K_x,dim_image*dim_image))
	eigen_vector[1,] = c(rep(c(rep(1,dim_image/3),rep(0,dim_image*2/3)),dim_image/3),rep(0,(dim_image*2/3)*dim_image))
	eigen_vector[2,] = c(rep(0,(dim_image*1/3)*dim_image),rep(c(rep(0,dim_image/3),rep(1,dim_image/3),rep(0,dim_image/3)),dim_image/3),rep(0,(dim_image*1/3)*dim_image))
	eigen_vector[3,] = c(rep(0,(dim_image*2/3)*dim_image),rep(c(rep(0,dim_image*2/3),rep(1,dim_image*1/3)),dim_image*1/3))

	####standarize the eigenimage
	for(jj in 1:K_x){
		eigen_vector[jj,] = eigen_vector[jj,]/sqrt(sum(eigen_vector[jj,]^2))
	}

	####we generate the eigen_score for each subject based on the normal distribution, following the mixed effect representation
	eigen_score = array(0,c(N_subj,K_x))

	####This is the variance for the noraml distributions generating the eigen_scores
	eigen_var = rep(0,K_x)
	eigen_var[1:K_x] = 0.5^((1:K_x) - 1)
	eigen_sd = sqrt(eigen_var)
	for(jj in 1:K_x){
		eigen_score[,jj] = rnorm(N_subj, mean = 0, sd = eigen_sd[jj])
	}

	####generate the images, X
	true.funcs = eigen_score%*%eigen_vector

	####demean the process, get the mu(t) out of the X
	mean.funcs = apply(true.funcs,2,mean)

	for(ii in 1:N_subj){
		true.funcs[ii,] = true.funcs[ii,] - mean.funcs
	}

	####fast svd to get eigenvalues and eigenvectors 
	eigen_svd = fast.svd(t(true.funcs),tol = 0.0001)

	####estimated eigenimages
	eigenimage_est = t(eigen_svd$u)

	####whether to keep the direction of the eigenimages
	ind_rev = rep(0,K_x)

	####the indicator of whehter we need to reverse the eigenimage to fit the direction of our simulated data
	for(jj in 1:K_x){
	   if(sum(eigenimage_est[jj,]*eigen_vector[jj,])<0){
		  eigenimage_est[jj,] = - eigenimage_est[jj,]
			ind_rev[jj] = 1
	  }

	}
	
	##### define the true beta function
	true_Beta = beta_true[1]*eigen_vector[1,] + beta_true[2]*eigen_vector[2,] + beta_true[3]*eigen_vector[3,] 

	##### define missing mechanism function
	true_Beta_R = beta_r_true[1]*eigen_vector[1,] + beta_r_true[2]*eigen_vector[2,] + beta_r_true[3]*eigen_vector[3,]

	#####generate covariates in the functional regressions and missing parts
	X_cov = array(0,dim = c(N_subj, q_cov))
	X_cov[,1] = runif(N_subj)
	X_cov[,2] = rnorm(N_subj, mean = 0, sd = 1)
	X_cov[,3] = rbinom(N_subj,1,0.5)
	X_cov[X_cov[,3]==0,3] = -1

	#####generate fully observed outcomes	
	outcomes <- sapply(1:N_subj, function(u) sum(true.funcs[u,]*true_Beta))+rnorm(N_subj, 0, sqrt(sigma2_delta_true)) + gamma_true[1]*X_cov[,1] + gamma_true[2]*X_cov[,2] + gamma_true[3]*X_cov[,3] + alpha_true

	##generate the missing data
	temp_unif = runif(N_subj,0,1)
	temp_prod = sapply(1:N_subj, function(u) sum(true.funcs[u,]*true_Beta_R)) + phi_true*outcomes  + gamma_r_true[1]*X_cov[,1] + gamma_r_true[2]*X_cov[,2] + alpha_r_true
	fun_mis<-function(x) exp(x)/(1 + exp(x))
	temp_prod = fun_mis(temp_prod)

	#####non-ignorable non-response
	R = rbinom(N_subj,1,temp_prod)
	xi = t(t(eigen_svd$v)*eigen_svd$d)
	Y = rep(NA,N_subj)
	Y[R==0]=outcomes[R==0]
		
	#####write the data to files
	if(CIR == 1){
		write(outcomes,file=Outcome_file,ncol=1,append=F)
		write(Y,file=Y_file,ncol=1,append=F)
		write(X_cov, file=Cov_file,ncol=q_cov,append=F)
		write(R,file = R_file,ncol=1,append=F)
		write(xi,file = Xi_file,ncol=K_x,append=F)
		write(ind_rev,file = Eigen_ind,ncol=K_x,append=F)
	}else{	
		write(outcomes,file=Outcome_file,ncol=1,append=T)
		write(Y,file=Y_file,ncol=1,append=T)
		write(X_cov, file=Cov_file,ncol=q_cov,append=T)
		write(R,file = R_file,ncol=1,append=T)
		write(xi,file = Xi_file,ncol=K_x,append=T)
		write(ind_rev,file = Eigen_ind,ncol=K_x,append=T)
	}
}


####analyze the simulated data
for(CIR in 1:N_REPS){
	####read data
	Y = matrix(scan(Y_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1)   ##outcome subject to missingness
	outcomes = matrix(scan(Outcome_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1)  ##true outcomes
	X_cov = matrix(scan(Cov_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=q_cov)  ##covariates in SoI regression
	R = matrix(scan(R_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1)  ##missing indicator
	xi = matrix(scan(Xi_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=K_x)  ##eigenscores from FPCA, and the FPCA could be easily conducted through matlab
	ind_rev = scan(Eigen_ind,skip=(CIR-1)*1,nlines=1)  ##for correcting the sign of the estimation in the simulation studies
	cat("\n")
	print(paste("Replication:", CIR))
	
	print("#########True values")
	print("####BSOI model:")
	print("##beta:")
	print(beta_true)
	print("##gamma:")
	print(gamma_true)
	print("####Tilting model:")
	print("##beta_r:")
	print(beta_r_true)
	print("##gamma_r:")
	print(gamma_r_true)
	print("##phi:")
	print(phi_true)
	
	#covariates in the tilting model
	X_cov_NN = X_cov[,ind_cov_NN]
	N_mis = sum(R==1)
	
	####SoI_IN procedure
	cat("\n")
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
	model_IN1 <- BSOIIN(Y, R, xi, X_cov, n_iter, alpha_0, sigma2_alpha0,  beta_k0, sigma2_betak0,   gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set1, beta_set1,  gamma_cov_set1, sigma2_delta_set1)
	model_IN2 <- BSOIIN(Y, R, xi, X_cov, n_iter, alpha_0, sigma2_alpha0,  beta_k0, sigma2_betak0,   gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set2, beta_set2,  gamma_cov_set2, sigma2_delta_set2)
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
	
    ##record results
	Ealpha_IN[CIR,] = mean(c(model_IN1$alpha[(n_burnin+1):n_iter],model_IN2$alpha[(n_burnin+1):n_iter],model_IN3$alpha[(n_burnin+1):n_iter]))
	Ebeta_IN[CIR,] = apply(rbind(model_IN1$beta[(n_burnin+1):n_iter,],model_IN2$beta[(n_burnin+1):n_iter,], model_IN3$beta[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
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
	var_ym =20
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
	Ebeta_NN[CIR,] = apply(rbind(model_NN1$beta[(n_burnin+1):n_iter,],model_NN2$beta[(n_burnin+1):n_iter,],model_NN3$beta[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
	Egamma_NN[CIR,] = apply(rbind(model_NN1$gamma_cov[(n_burnin+1):n_iter,],model_NN2$gamma_cov[(n_burnin+1):n_iter,],model_NN3$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_NN[CIR,] = mean(c(model_NN1$sigma2_delta[(n_burnin+1):n_iter],model_NN2$sigma2_delta[(n_burnin+1):n_iter],model_NN3$sigma2_delta[(n_burnin+1):n_iter]))
	
	Ealpha_r_NN[CIR,] = mean(c(model_NN1$alpha_r[(n_burnin+1):n_iter],model_NN2$alpha_r[(n_burnin+1):n_iter],model_NN3$alpha_r[(n_burnin+1):n_iter]))
	Ebeta_r_NN[CIR,] = apply(rbind(model_NN1$beta_r[(n_burnin+1):n_iter,],model_NN2$beta_r[(n_burnin+1):n_iter,],model_NN3$beta_r[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
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
	
	##print the results
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
	
	####SoI_Full procedure
	cat("\n")
	print("BSOI_Full:")
	Y = outcomes  ##This procedure assumes the response is fully observed
	
	##hyperparameters
	alpha_0 = 0
	sigma2_alpha0 = 100
	beta_k0 = rep(0,K_x)
	sigma2_betak0 = rep(100,K_x)
	gamma_q0 = rep(0.01,q_cov)
	sigma2_gammaq0 = rep(100,q_cov)
	a_delta0 = 9
	b_delta0 = 3
    
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
    
    ##conduct the Bayesian analysis using the function SoIFull()
    ##Only the data inputs are compulsory, and we may implement the analysis as
    ##model_Full <- BSOIFull(Y, xi, X_cov, n_iter)	
	model_Full1 <- BSOIFull(Y, xi, X_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0,  sigma2_betak0,  gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set1, beta_set1, gamma_cov_set1, sigma2_delta_set1)
	model_Full2 <- BSOIFull(Y, xi, X_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0,  sigma2_betak0,  gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set2, beta_set2, gamma_cov_set2, sigma2_delta_set2)
	model_Full3 <- BSOIFull(Y, xi, X_cov, n_iter, alpha_0, sigma2_alpha0, beta_k0,  sigma2_betak0,  gamma_q0, sigma2_gammaq0, a_delta0, b_delta0,  alpha_set3, beta_set3, gamma_cov_set3, sigma2_delta_set3)
	
	#beta and gamma are major parameters, and the convergence is checked through the traceplots of multiple chains
	par(mfrow=c(3,2))
	for(iii in 1:K_x){
		plot(model_Full1$beta[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("beta",iii), main = "BSOIFull" )
		lines(model_Full2$beta[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
		lines(model_Full3$beta[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
	par(mfrow=c(4,2))
	for(iii in 1:q_cov){
	  plot(model_Full1$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "red", xlab = "Iterations after burn-in", ylab = paste("gamma",iii), main = "BSOIFull" )
	  lines(model_Full2$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "blue")
	  lines(model_Full3$gamma_cov[(n_burnin+1):n_iter,iii],type = "l",col = "green")
	}
	
    ##record results
	Ealpha_Full[CIR,] = mean(c(model_Full1$alpha[(n_burnin+1):n_iter],model_Full2$alpha[(n_burnin+1):n_iter],model_Full3$alpha[(n_burnin+1):n_iter]))
	Ebeta_Full[CIR,] = apply(rbind(model_Full1$beta[(n_burnin+1):n_iter,],model_Full2$beta[(n_burnin+1):n_iter,], model_Full3$beta[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
	Egamma_Full[CIR,] = apply(rbind(model_Full1$gamma_cov[(n_burnin+1):n_iter,],model_Full2$gamma_cov[(n_burnin+1):n_iter,],model_Full3$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_Full[CIR,] = mean(c(model_Full1$sigma2_delta[(n_burnin+1):n_iter],model_Full2$sigma2_delta[(n_burnin+1):n_iter],model_Full3$sigma2_delta[(n_burnin+1):n_iter]))
	
	SDalpha_Full[CIR,] = sd(c(model_Full1$alpha[(n_burnin+1):n_iter],model_Full2$alpha[(n_burnin+1):n_iter],model_Full3$alpha[(n_burnin+1):n_iter]))
	SDbeta_Full[CIR,] = apply(rbind(model_Full1$beta[(n_burnin+1):n_iter,],model_Full2$beta[(n_burnin+1):n_iter,], model_Full3$beta[(n_burnin+1):n_iter,]),2,sd)
	SDgamma_Full[CIR,] = apply(rbind(model_Full1$gamma_cov[(n_burnin+1):n_iter,],model_Full2$gamma_cov[(n_burnin+1):n_iter,],model_Full3$gamma_cov[(n_burnin+1):n_iter,]),2,sd)
	SDdelta_Full[CIR,] = sd(c(model_Full1$sigma2_delta[(n_burnin+1):n_iter],model_Full2$sigma2_delta[(n_burnin+1):n_iter],model_Full3$sigma2_delta[(n_burnin+1):n_iter]))
	
	##print results
	cat("\n")
	print("#########Estimation results of BSOIFull:")
	print("##beta:")
	print("mean:")
	print(Ebeta_Full[CIR,])
	print("SD:")
	print(SDbeta_Full[CIR,])
	print("##gamma:")
	print("mean:")
	print(Egamma_Full[CIR,])
	print("SD:")
	print(SDgamma_Full[CIR,])
} 