set.seed(as.numeric(Sys.time()))
library(BSOINN)  ##load the library

####file names
Y_file = "Y.txt"
Cov_file = "COV.txt"  
R_file = "R.txt"
Xi_file = "XI.txt"

####define constants
N_REPS = 5
N_subj=802     ##sample size
n_iter=10000    ##total iterations
n_burnin=5000  ##burn-in phase
K_x = 5        ##number of eigenimages retained   
q_cov = 7      ##dimension of covariates in SoI
q_cov_NN = 6   ##dimension of covariates in missing mechanism after removing the instrument 
ind_cov_NN = c(T,T,F,T,T,T,T)  ##we define educational level (the third) as an instrument


###define objects for prediction results
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

R = matrix(scan(R_file,skip=0,nlines=N_subj),nrow=N_subj,ncol=1,byrow = T)  ##missing indicator

####gendata
for(CIR in 1:N_REPS){
	##Prediction
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


####read the values
Y = matrix(scan(Y_file,skip=0,nlines=N_subj),nrow=N_subj,ncol=1, byrow = T)   ##outcome subject to missingness
X_cov = matrix(scan(Cov_file,skip=0,nlines=N_subj),nrow=N_subj,ncol=q_cov, byrow = T)  ##covariates in SoI regression
R = matrix(scan(R_file,skip=0,nlines=N_subj),nrow=N_subj,ncol=1,byrow = T)  ##missing indicator
xi = matrix(scan(Xi_file,skip=0,nlines=N_subj),nrow=N_subj,ncol=K_x, byrow=T)  ##eigenscores from FPCA, and the FPCA could be easily conducted through matlab or the fast.svd function in R
	
X_cov_NN = X_cov[,ind_cov_NN]
N_mis = sum(R==1)
posi_mis = which(R==1)
posi_obs = which(R==0)
Y_true_file = "Y_true.txt"
outcomes = matrix(scan(Y_true_file,skip=0,nlines=N_subj),nrow=N_subj,ncol=1, byrow = T)

for(CIR in 1:N_REPS){
	cat("\n")
	print(paste("Replication:", CIR))

	##Prediction
	##Step 1, separate the data set into training set and test set
	posi_train = scan("posi_train.txt",skip=N_subj/2*(CIR-1),nlines=N_subj/2)
	posi_test = scan("posi_test.txt",skip=N_subj/2*(CIR-1),nlines=N_subj/2)

	Y_train = Y[posi_train]
	R_train = R[posi_train]
	mean_Y_train = mean(Y_train[R_train == 0])
	sd_Y_train = sd(Y_train[R_train == 0])
	Y_train[R_train == 0] =  (Y_train[R_train == 0] - mean_Y_train)/sd_Y_train 
	xi_train  = xi[posi_train,]
	X_cov_train = X_cov[posi_train,]
	X_cov_NN_train = X_cov_NN[posi_train,]
	train_n_mis = sum(R_train == 1)
	
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
	Ebeta_IN_train[CIR,] = apply(rbind(model_IN1_train$beta[(n_burnin+1):n_iter,],model_IN2_train$beta[(n_burnin+1):n_iter,], model_IN3_train$beta[(n_burnin+1):n_iter,]),2,mean)
	Egamma_IN_train[CIR,] = apply(rbind(model_IN1_train$gamma_cov[(n_burnin+1):n_iter,],model_IN2_train$gamma_cov[(n_burnin+1):n_iter,],model_IN3_train$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_IN_train[CIR,] = mean(c(model_IN1_train$sigma2_delta[(n_burnin+1):n_iter],model_IN2_train$sigma2_delta[(n_burnin+1):n_iter],model_IN3_train$sigma2_delta[(n_burnin+1):n_iter]))
	
	##test set prediction	
	Y_test_true = outcomes[posi_test]

	R_test = R[posi_test]
	xi_test  = xi[posi_test,]
	X_cov_test = X_cov[posi_test,]
	
	Y_test_pred_IN = xi_test%*%as.matrix(Ebeta_IN_train[CIR,]) + X_cov_test%*%as.matrix(Egamma_IN_train[CIR,])

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
	Ebeta_NN_train[CIR,] = apply(rbind(model_NN1_train$beta[(n_burnin+1):n_iter,],model_NN2_train$beta[(n_burnin+1):n_iter,],model_NN3_train$beta[(n_burnin+1):n_iter,]),2,mean)
	Egamma_NN_train[CIR,] = apply(rbind(model_NN1_train$gamma_cov[(n_burnin+1):n_iter,],model_NN2_train$gamma_cov[(n_burnin+1):n_iter,],model_NN3_train$gamma_cov[(n_burnin+1):n_iter,]),2,mean)
	Edelta_NN_train[CIR,] = mean(c(model_NN1_train$sigma2_delta[(n_burnin+1):n_iter],model_NN2_train$sigma2_delta[(n_burnin+1):n_iter],model_NN3_train$sigma2_delta[(n_burnin+1):n_iter]))

	Ealpha_r_NN_train[CIR,] = mean(c(model_NN1_train$alpha_r[(n_burnin+1):n_iter],model_NN2_train$alpha_r[(n_burnin+1):n_iter],model_NN3_train$alpha_r[(n_burnin+1):n_iter]))
	Ebeta_r_NN_train[CIR,] = apply(rbind(model_NN1_train$beta_r[(n_burnin+1):n_iter,],model_NN2_train$beta_r[(n_burnin+1):n_iter,],model_NN3_train$beta_r[(n_burnin+1):n_iter,]),2,mean)
	Egamma_r_NN_train[CIR,] = apply(rbind(model_NN1_train$gamma_cov_r[(n_burnin+1):n_iter,],model_NN2_train$gamma_cov_r[(n_burnin+1):n_iter,],model_NN3_train$gamma_cov_r[(n_burnin+1):n_iter,]),2,mean)
	Ephi_NN_train[CIR,] = mean(c(model_NN1_train$phi[(n_burnin+1):n_iter],model_NN2_train$phi[(n_burnin+1):n_iter],model_NN3_train$phi[(n_burnin+1):n_iter]))

	Y_test_pred_NN =  xi_test%*%as.matrix(Ebeta_NN_train[CIR,]) + X_cov_test%*%as.matrix(Egamma_NN_train[CIR,])
	cat("\n")
	print("BSOI_NN Prediction:")
	Acc_pred_NN[CIR] = cor(Y_test_pred_NN,Y_test_true)
	print(Acc_pred_NN[CIR])
} 
