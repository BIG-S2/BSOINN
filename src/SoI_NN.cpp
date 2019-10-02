// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List BSOINN(vec Y, vec R, mat xi, mat Z_cov, mat Z_cov_NN, int n_iter, double alpha_0, double sigma2_alpha0, vec beta_k0, vec sigma2_betak0,  vec gamma_q0, vec sigma2_gammaq0,double a_delta0, double b_delta0, vec varphi_0, mat sigma_varphi0, double var_ym, double var_varphi, double alpha, vec beta, vec gamma_cov, double sigma2_delta, double alpha_r, vec beta_r, vec gamma_cov_r, double phi) {
	
	
	 // Initializing
	 int N_subj = Y.n_elem;
	 int N_mis = sum(R);
	 int K_x = xi.n_cols;
	 int q_cov = Z_cov.n_cols;
	 int q_cov_NM = Z_cov_NN.n_cols;
	 
	 
	 vec data_alpha(n_iter);
	 mat data_beta(n_iter,K_x);
	 mat data_gamma_cov(n_iter,q_cov);
	 vec data_sigma2_delta(n_iter);
	 vec data_alpha_r(n_iter);
	 mat data_beta_r(n_iter,K_x);
	 mat data_gamma_cov_r(n_iter,q_cov_NM);
	 vec data_phi(n_iter);
	 ivec acce_ym_n = zeros<ivec>(N_mis);
	 int acce_varphi_n = 0;
	 
	 double sigma_beta_temp = 0;
	 double mean_beta_temp = 0;
	 double sigma_gamma_temp = 0;
	 double mean_gamma_temp = 0;
	 vec xi_sum_square = zeros<vec>(K_x);
	 vec XCov_sum_square = zeros<vec>(q_cov);
	 vec sigma_ym_temp = zeros<vec>(N_mis);
	 vec xi_beta_r = zeros<vec>(N_mis);
	 vec xi_beta = zeros<vec>(N_mis);
	 vec ym_new = zeros<vec>(N_mis);
	 vec likelihood1_ym = zeros<vec>(N_mis);
	 vec likelihood2_ym = zeros<vec>(N_mis);
	 vec prob_ym = zeros<vec>(N_mis);
	 int current_posi = 0;
	 uvec xi_col_indices(K_x - 1);
	 uvec gamma_col_indices(q_cov - 1);
	 uvec missing_ind(N_mis);
	 
	 mat ee_temp = diagmat(ones<vec>(K_x+q_cov_NM+2));
	 mat ee_temp_kk = zeros<mat>(K_x+q_cov_NM+2,1);
	 vec varphi_new(K_x + q_cov_NM+2);
	 vec varphi_temp(K_x + q_cov_NM+2);
	 mat ee_temp_matrix = ones<mat>(N_subj,K_x+q_cov_NM+2);
	 
	 ee_temp_matrix(span(0,N_subj-1),span(1,q_cov_NM)) = Z_cov_NN;
     ee_temp_matrix(span(0,N_subj-1),span(1+q_cov_NM,q_cov_NM+K_x)) = xi;
	 
	 double likelihood1_varphi;
	 double likelihood2_varphi;
	 mat inv_sigma_varphi0 = inv_sympd(sigma_varphi0);
	 double prob_varphi;
	
	 for(int i = 0; i < N_subj; i++){
		 if(R[i] == 1){
			missing_ind[current_posi] = i;
			current_posi++;
		 }
	 }
	 
	 for (int kk =0; kk < K_x; kk++){
		xi_sum_square[kk] = sum(pow(norm(xi.col(kk)),2));
	    }
		
	 for (int kk =0; kk < q_cov; kk++){
		XCov_sum_square[kk] = sum(pow(norm(Z_cov.col(kk)),2));
	    }
		
	 varphi_temp[0] = alpha_r;
	 varphi_temp(span(1,q_cov_NM)) = gamma_cov_r;
     varphi_temp(span(q_cov_NM+1,q_cov_NM+K_x)) = beta_r;
	 varphi_temp[1+K_x+q_cov_NM] = phi;
		
	 
	 for (int it = 0; it < n_iter; it++){
		 
		 alpha = (sum(Y - xi*beta -Z_cov*gamma_cov)/sigma2_delta + alpha_0/sigma2_alpha0)/(1/sigma2_alpha0 + N_subj/sigma2_delta) + rnorm(1,0,1)[0]*sqrt(1/(1/sigma2_alpha0 + N_subj/sigma2_delta));
		 //alpha = sum(Y - xi*beta -Z_cov*gamma_cov)/N_subj + random_alpha[it]*sqrt(sigma2_delta/N_subj);
		 
		 for (int kk =0; kk < K_x; kk++){
			 current_posi = 0;
			 sigma_beta_temp = 1/(xi_sum_square[kk]/sigma2_delta + 1/sigma2_betak0[kk]);
			 if(K_x ==1){
				mean_beta_temp = sigma_beta_temp*(beta_k0[kk]/sigma2_betak0[kk] + sum((Y - alpha - Z_cov*gamma_cov)%xi.col(kk))/sigma2_delta);
			 }else{
				 for(int kk2 = 0; kk2 < K_x; kk2++){
					 if(kk2!=kk){
						xi_col_indices[current_posi] = kk2;
						current_posi++;
					 }
				 }
				 mean_beta_temp = sigma_beta_temp*(beta_k0[kk]/sigma2_betak0[kk] + sum((Y - alpha - Z_cov*gamma_cov - xi.cols(xi_col_indices)*beta(xi_col_indices))%xi.col(kk))/sigma2_delta);
			 }
			 beta[kk] = mean_beta_temp + sqrt(sigma_beta_temp)*rnorm(1,0,1)[0];
			 //beta[kk] = mean_beta_temp + sqrt(sigma_beta_temp)*random_beta(kk,it);
		 }
		 
		 for (int kk = 0; kk < q_cov; kk++){
			 current_posi = 0;
			 sigma_gamma_temp = 1/(XCov_sum_square[kk]/sigma2_delta + 1/sigma2_gammaq0[kk]);
			  if(q_cov ==1){
				mean_gamma_temp = sigma_gamma_temp*(gamma_q0[kk]/sigma2_gammaq0[kk] + sum((Y - alpha - xi*beta)%Z_cov.col(kk))/sigma2_delta);
			 }else{
				 for(int kk2 = 0; kk2 < q_cov; kk2++){
					 if(kk2!=kk){
						gamma_col_indices[current_posi] = kk2;
						current_posi++;
					 }
				 }
				 mean_gamma_temp = sigma_gamma_temp*(gamma_q0[kk]/sigma2_gammaq0[kk] + sum((Y - alpha - xi*beta - Z_cov.cols(gamma_col_indices)*gamma_cov(gamma_col_indices))%Z_cov.col(kk))/sigma2_delta);
			 }
			 gamma_cov[kk] = mean_gamma_temp + sqrt(sigma_gamma_temp)*rnorm(1,0,1)[0];
			 //gamma_cov[kk] = mean_gamma_temp + sqrt(sigma_gamma_temp)*random_gamma(kk,it);
		 }
		 
		 
		 sigma2_delta = 1 / rgamma(1, a_delta0 + 0.5 * N_subj, 1 / (b_delta0 + 0.5*sum(pow(norm(Y - alpha - xi*beta -Z_cov*gamma_cov),2))))[0];
		 
		 //update y_mis
		
		 xi_beta_r = alpha_r + Z_cov_NN.rows(missing_ind)*gamma_cov_r + xi.rows(missing_ind)*beta_r;
		 
		 xi_beta = alpha + Z_cov.rows(missing_ind)*gamma_cov + xi.rows(missing_ind)*beta;
		 
		 sigma_ym_temp = sqrt(var_ym/(1/sigma2_delta + pow(phi,2)*exp(xi_beta_r)/pow(exp(xi_beta_r) + 1,2)));
		 
		 ym_new = as<arma::vec>(rnorm(N_mis))%sigma_ym_temp + Y(missing_ind);
		 //ym_new = random_y(span(0,N_mis-1),it)%sigma_ym_temp + Y(missing_ind);
		 
		 likelihood1_ym = pow(Y(missing_ind) - xi_beta,2)/sigma2_delta - 2*phi*Y(missing_ind) + 2*log(1 + exp(xi_beta_r + phi*Y(missing_ind)));
		 
		 likelihood2_ym = pow(ym_new - xi_beta,2)/sigma2_delta - 2*phi*ym_new + 2*log(1 + exp(xi_beta_r + phi*ym_new));
		 
		 prob_ym = exp(-0.5*(likelihood2_ym - likelihood1_ym));
		 
		 for(int kk =0; kk < N_mis; kk++){
			 if(runif(1)[0]<prob_ym[kk]){
			   //if(random_u_y(kk,it)<prob_ym[kk]){ 
				 Y[missing_ind[kk]] = ym_new[kk];
				 acce_ym_n[kk]++;
			 }
		 }
		
		 //update varphi, remember the elements of matrix should be ().
		 
		 ee_temp = inv_sigma_varphi0;
	     
		 for(int kk =0; kk <N_subj; kk++){
			 ee_temp_kk(0,0) = 1;
			 ee_temp_kk(span(1,q_cov_NM),0) = Z_cov_NN.row(kk).t();
			 ee_temp_kk(span(q_cov_NM+1,q_cov_NM+K_x),0) = xi.row(kk).t();
			 ee_temp_kk(q_cov_NM+K_x+1,0) = Y[kk] ;
			 
			 ee_temp = ee_temp + 0.25*(ee_temp_kk*(ee_temp_kk.t()));
		 }
		 
		 ee_temp = inv_sympd(ee_temp);
		 ee_temp = var_varphi*ee_temp;
		
		
		 varphi_new = (chol(ee_temp).t())*(as<arma::vec>(rnorm(2+q_cov_NM+K_x))) + varphi_temp;
		 
		 //varphi_new = chol(ee_temp).t()*random_var(span(0,6),it) + varphi_temp; //varphi_temp;//chol(ee_temp)*random_number(span(0,6),it) + varphi_temp; //mvrnorm(1, varphi_temp, ee_temp).row(0).t(); //chol(ee_temp)*(as<arma::vec>(rnorm(2+q_cov_NM+K_x))) + varphi_temp; 
		
		 		 
		 ee_temp_matrix(span(0,N_subj-1),1+q_cov_NM+K_x) = Y;
		 
		 likelihood1_varphi = -2*sum(R%(ee_temp_matrix*varphi_temp)) + sum(0.5*((varphi_temp - varphi_0).t())*inv_sigma_varphi0*(varphi_temp - varphi_0)) + 2*sum(log(1 + exp(ee_temp_matrix*varphi_temp)));
		 
		 likelihood2_varphi = -2*sum(R%(ee_temp_matrix*varphi_new)) + sum(0.5*((varphi_new - varphi_0).t())*inv_sigma_varphi0*(varphi_new - varphi_0)) + 2*sum(log(1 + exp(ee_temp_matrix*varphi_new)));
		 
		 prob_varphi = exp(-0.5*(likelihood2_varphi - likelihood1_varphi));
		 
		 if(runif(1)[0]<prob_varphi){
		   //if(random_u_var[it]<prob_varphi){
			 alpha_r = varphi_new[0];
			 gamma_cov_r = varphi_new(span(1,q_cov_NM));
			 beta_r = varphi_new(span(q_cov_NM+1,q_cov_NM+K_x));
			 phi = varphi_new[1+q_cov_NM+K_x];
			 acce_varphi_n++;
			 varphi_temp = varphi_new;
		 }
		 	 	
		data_alpha[it] = alpha;
		data_beta.row(it) = beta.t();
		data_gamma_cov.row(it) = gamma_cov.t();
		data_sigma2_delta[it] = sigma2_delta;
		data_alpha_r[it] = alpha_r;
		data_beta_r.row(it) = beta_r.t();
		data_gamma_cov_r.row(it) = gamma_cov_r.t();
		data_phi[it] = phi;
	
 
		if (it > 0 && it % (n_iter/100) == 0)
			Rprintf("\rRunning MCMC, %d%% completed...", it / (n_iter/100));
	 }		  
			  
	 return List::create(Named("alpha") = data_alpha,
                         Named("beta") = data_beta,
                         Named("sigma2_delta") = data_sigma2_delta,
                         Named("gamma_cov") = data_gamma_cov,
						 Named("alpha_r") = data_alpha_r,
                         Named("beta_r") = data_beta_r,
                         Named("phi") = data_phi,
                         Named("gamma_cov_r") = data_gamma_cov_r,
						 Named("acce_ym") = acce_ym_n,
						 Named("acce_varphi_n") = acce_varphi_n);

//Please note that the elements of a matrix is not [] but (), and chol() in arma -> trans(chol()) in R				 

}

