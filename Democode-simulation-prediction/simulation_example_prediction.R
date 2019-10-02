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
N_subj=500     ##sample size
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
Ealpha_IN_train = array(0, dim = c(N_REPS,1))
Egamma_IN_train = array(0,dim = c(N_REPS,q_cov))
Ebeta_IN_train = array(0, dim = c(N_REPS,K_x))
Edelta_IN_train = array(0, dim = c(N_REPS,1))

Ealpha_NN_train = array(0, dim = c(N_REPS,1))
Egamma_NN_train = array(0,dim = c(N_REPS,q_cov))
Ebeta_NN_train = array(0, dim = c(N_REPS,K_x))
Edelta_NN_train = array(0, dim = c(N_REPS,1))
Ealpha_r_NN_train = array(0, dim = c(N_REPS,1))
Egamma_r_NN_train = array(0,dim = c(N_REPS,q_cov_NN))
Ebeta_r_NN_train =array(0, dim = c(N_REPS,K_x))
Ephi_NN_train = array(0, dim = c(N_REPS,1))

Ealpha_Full_train = array(0, dim = c(N_REPS,1))
Egamma_Full_train = array(0,dim = c(N_REPS,q_cov))
Ebeta_Full_train = array(0, dim = c(N_REPS,K_x))
Edelta_Full_train = array(0, dim = c(N_REPS,1))

Acc_pred_IN = rep(0,N_REPS)
Acc_pred_NN = rep(0,N_REPS)
Acc_pred_Full = rep(0,N_REPS)

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
	Y = rep(999,N_subj)
	Y[R==0]=outcomes[R==0]
		
	#####write the data to files
	if(CIR == 1){
		write(outcomes,file=Outcome_file,ncol=1,append=F)
		write(Y,file=Y_file,ncol=1,append=F)
		write(t(X_cov), file=Cov_file,ncol=q_cov,append=F)
		write(R,file = R_file,ncol=1,append=F)
		write(t(xi),file = Xi_file,ncol=K_x,append=F)
		write(ind_rev,file = Eigen_ind,ncol=K_x,append=F)
	}else{	
		write(outcomes,file=Outcome_file,ncol=1,append=T)
		write(Y,file=Y_file,ncol=1,append=T)
		write(t(X_cov), file=Cov_file,ncol=q_cov,append=T)
		write(R,file = R_file,ncol=1,append=T)
		write(t(xi),file = Xi_file,ncol=K_x,append=T)
		write(ind_rev,file = Eigen_ind,ncol=K_x,append=T)
	}
	
	##Step 1, separate the data set into training set and test set
posi_mis = which(R==1)
posi_obs = which(R==0)

train_n_mis = floor(sum(R==1)/2)
posi_train_mis_posi = sample(x = sum(R==1),size = train_n_mis, rep = F)
posi_train_mis = posi_mis[posi_train_mis_posi]
train_n_obs = N_subj/2 - train_n_mis
posi_train_obs_posi = sample(x = sum(R==0),size = train_n_obs, rep = F)
posi_train_obs = posi_obs[posi_train_obs_posi]
posi_train = c(posi_train_mis, posi_train_obs)

test_n_mis = sum(R==1) - train_n_mis
posi_test_mis = posi_mis[-posi_train_mis_posi]
test_n_obs = N_subj/2 - test_n_mis
posi_test_obs = posi_obs[-posi_train_obs_posi]
posi_test = c(posi_test_mis,posi_test_obs)

	if(CIR == 1){
		write(posi_train,file="posi_train.txt",ncol=1,append=F)
		write(posi_test,file="posi_test.txt",ncol=1,append=F)	
	}else{
		write(posi_train,file="posi_train.txt",ncol=1,append=T)
		write(posi_test,file="posi_test.txt",ncol=1,append=T)
	}
	
}

####analyze the simulated data
for(CIR in 1:N_REPS){
	####read data
	Y = matrix(scan(Y_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1)   ##outcome subject to missingness
	outcomes = matrix(scan(Outcome_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1)  ##true outcomes
	X_cov = matrix(scan(Cov_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=q_cov,byrow = T)  ##covariates in SoI regression
	R = matrix(scan(R_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=1)  ##missing indicator
	xi = matrix(scan(Xi_file,skip=(CIR-1)*N_subj,nlines=N_subj),nrow=N_subj,ncol=K_x, byrow = T)  ##eigenscores from FPCA, and the FPCA could be easily conducted through matlab
	ind_rev = scan(Eigen_ind,skip=(CIR-1)*1,nlines=1)  ##for correcting the sign of the estimation in the simulation studies
	cat("\n")
	print(paste("Replication:", CIR))
	
	
	posi_train = scan("posi_train.txt",skip=(CIR-1)*1,nlines=N_subj/2)
	posi_test = scan("posi_test.txt",skip=(CIR-1)*1,nlines=N_subj/2)

	Y_train = Y[posi_train]
	R_train = R[posi_train]
	xi_train  = xi[posi_train,]
	X_cov_train = X_cov[posi_train,]
	X_cov_NN_train = X_cov_train[,ind_cov_NN]

	##test set prediction	
	Y_train_true = outcomes[posi_train]

	Y_test_true = outcomes[posi_test]
	R_test = R[posi_test]
	xi_test  = xi[posi_test,]
	X_cov_test = X_cov[posi_test,]
	train_n_mis = sum(R_train==1)
		
	##Prediction
	##training set evaluation
	##temp values for missing components
	Y_train[R_train==1] = rnorm(train_n_mis,mean = 0, sd = 1)
	
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
	model_IN1_train <- BSOIIN(Y_train, R_train, xi_train, X_cov_train, n_iter)
	model_IN2_train <- BSOIIN(Y_train, R_train, xi_train, X_cov_train, n_iter)
	model_IN3_train <- BSOIIN(Y_train, R_train, xi_train, X_cov_train, n_iter)
	
	##record results
	Ealpha_IN_train[CIR,] = mean(c(model_IN1_train$alpha[(n_burnin+1):n_iter],model_IN2_train$alpha[(n_burnin+1):n_iter],model_IN3_train$alpha[(n_burnin+1):n_iter]))
	Ebeta_IN_train[CIR,] = apply(rbind(model_IN1_train$beta[(n_burnin+1):n_iter,],model_IN2_train$beta[(n_burnin+1):n_iter,], model_IN3_train$beta[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
	Egamma_IN_train[CIR,] = apply(rbind(model_IN1_train$gamma_cov[(n_burnin+1):n_iter,],model_IN2_train$gamma_cov[(n_burnin+1):n_iter,],model_IN3_train$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_IN_train[CIR,] = mean(c(model_IN1_train$sigma2_delta[(n_burnin+1):n_iter],model_IN2_train$sigma2_delta[(n_burnin+1):n_iter],model_IN3_train$sigma2_delta[(n_burnin+1):n_iter]))



	Y_test_pred_IN =  xi_test%*%(as.matrix(Ebeta_IN_train[CIR,])*(1 - 2*ind_rev)) + X_cov_test%*%as.matrix(Egamma_IN_train[CIR,])

	cat("\n")
	print("BSOI_IN Prediction:")
	Acc_pred_IN[CIR] = cor(Y_test_pred_IN,Y_test_true)
	print(Acc_pred_IN[CIR])

	####SoI_NN procedure
	##temp values for missing components
	Y_train[R_train==1] = rnorm(train_n_mis,mean = 0, sd = 1)

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
	model_NN1_train <- BSOINN(Y_train, R_train, xi_train, X_cov_train, X_cov_NN_train, n_iter)
	model_NN2_train <- BSOINN(Y_train, R_train, xi_train, X_cov_train, X_cov_NN_train, n_iter)
	model_NN3_train <- BSOINN(Y_train, R_train, xi_train, X_cov_train, X_cov_NN_train, n_iter)

	##record results
	Ealpha_NN_train[CIR,] = mean(c(model_NN1_train$alpha[(n_burnin+1):n_iter],model_NN2_train$alpha[(n_burnin+1):n_iter],model_NN3_train$alpha[(n_burnin+1):n_iter]))
	Ebeta_NN_train[CIR,] = apply(rbind(model_NN1_train$beta[(n_burnin+1):n_iter,],model_NN2_train$beta[(n_burnin+1):n_iter,],model_NN3_train$beta[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
	Egamma_NN_train[CIR,] = apply(rbind(model_NN1_train$gamma_cov[(n_burnin+1):n_iter,],model_NN2_train$gamma_cov[(n_burnin+1):n_iter,],model_NN3_train$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_NN_train[CIR,] = mean(c(model_NN1_train$sigma2_delta[(n_burnin+1):n_iter],model_NN2_train$sigma2_delta[(n_burnin+1):n_iter],model_NN3_train$sigma2_delta[(n_burnin+1):n_iter]))

	Ealpha_r_NN_train[CIR,] = mean(c(model_NN1_train$alpha_r[(n_burnin+1):n_iter],model_NN2_train$alpha_r[(n_burnin+1):n_iter],model_NN3_train$alpha_r[(n_burnin+1):n_iter]))
	Ebeta_r_NN_train[CIR,] = apply(rbind(model_NN1_train$beta_r[(n_burnin+1):n_iter,],model_NN2_train$beta_r[(n_burnin+1):n_iter,],model_NN3_train$beta_r[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
	Egamma_r_NN_train[CIR,] = apply(rbind(model_NN1_train$gamma_cov_r[(n_burnin+1):n_iter,],model_NN2_train$gamma_cov_r[(n_burnin+1):n_iter,],model_NN3_train$gamma_cov_r[(n_burnin+1):n_iter,]),2,mean)
	Ephi_NN_train[CIR,] = mean(c(model_NN1_train$phi[(n_burnin+1):n_iter],model_NN2_train$phi[(n_burnin+1):n_iter],model_NN3_train$phi[(n_burnin+1):n_iter]))

	Y_test_pred_NN =  xi_test%*%(as.matrix(Ebeta_NN_train[CIR,])*(1 - 2*ind_rev)) + X_cov_test%*%as.matrix(Egamma_NN_train[CIR,])
	cat("\n")
	print("BSOI_NN Prediction:")
	Acc_pred_NN[CIR] = cor(Y_test_pred_NN,Y_test_true)	
	print(Acc_pred_NN[CIR])
	
	###BSOI full prediction
	####SoI_Full procedure
	cat("\n")
	print("BSOI_Full:")
	Y_train = outcomes[posi_train]  ##This procedure assumes the response is fully observed
    
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
	model_Full1_train <- BSOIFull(Y_train, xi_train, X_cov_train, n_iter)
	model_Full2_train <- BSOIFull(Y_train, xi_train, X_cov_train, n_iter)
	model_Full3_train <- BSOIFull(Y_train, xi_train, X_cov_train, n_iter)

	
    ##record results
	Ealpha_Full_train[CIR,] = mean(c(model_Full1_train$alpha[(n_burnin+1):n_iter],model_Full2_train$alpha[(n_burnin+1):n_iter],model_Full3_train$alpha[(n_burnin+1):n_iter]))
	Ebeta_Full_train[CIR,] = apply(rbind(model_Full1_train$beta[(n_burnin+1):n_iter,],model_Full2_train$beta[(n_burnin+1):n_iter,], model_Full3_train$beta[(n_burnin+1):n_iter,]),2,mean)*(1 - 2*ind_rev)
	Egamma_Full_train[CIR,] = apply(rbind(model_Full1_train$gamma_cov[(n_burnin+1):n_iter,],model_Full2_train$gamma_cov[(n_burnin+1):n_iter,],model_Full3_train$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_Full_train[CIR,] = mean(c(model_Full1_train$sigma2_delta[(n_burnin+1):n_iter],model_Full2_train$sigma2_delta[(n_burnin+1):n_iter],model_Full3_train$sigma2_delta[(n_burnin+1):n_iter]))
	

	Y_test_pred_Full =  xi_test%*%(as.matrix(Ebeta_Full_train[CIR,]) *(1 - 2*ind_rev)) + X_cov_test%*%as.matrix(Egamma_Full_train[CIR,])
	
	cat("\n")
	print("BSOI_Full Prediction:")
	Acc_pred_Full[CIR] = cor(Y_test_pred_Full,Y_test_true)
	print(Acc_pred_Full[CIR])
} 