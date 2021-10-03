
functions{
// the Cholesky decomposition of AR1 matrix (analytical result)
 matrix AR1_L(int dim, real rho){
    matrix[dim, dim] ar1_mat_L = rep_matrix(0, dim, dim);
	real prev;
	
    ar1_mat_L[1,1] = 1;
    for(j in 2:dim)
    {
        ar1_mat_L[j,j] = sqrt(1 - rho^2);
    }
    
    for(j in 1:(dim-1))
    {
        for(i in (j+1):dim)
        {
            prev = ar1_mat_L[i-1,j];
			ar1_mat_L[i,j] = rho*prev;
        }
    }
    return(ar1_mat_L);
 }
}

data {
  int<lower=1> Ni;   // number of spatial units
  int<lower=1> Nk;   // number of groups
  int<lower=1> Nt;   // number of time point
  
  matrix[Ni, Ni] D;   //diagonal matrix with elements mi, where mi denotes the number of neighbors for the ith spatial unit
  matrix[Ni, Ni] W;   //adjacence matrix with elements 1 for Wij if i is a neighbor of j, 0 otherwise
  
  matrix[Ni, Nt] N_it;   // population size
  int Y[Ni, Nt];      // number of all cases
  int Y_severe[Ni, Nt];  // number of severe cases
  
  int Z_mild[Ni, Nk, Nt];  //number of mild cases for each serotype
  int Z_severe[Ni, Nk, Nt];   //number of severe cases for each serotype
  
  matrix[Ni, Nt] p_test_severe;    // empirical proportions of severe cases tested
  matrix[Ni, Nt] p_test_mild;   // empirical proportions of mild cases tested
  
  matrix[Nk, Nk] G_t_L;
  vector<lower = 0>[Nk] sigma_k; //standard deviation of cross-pathogen correlation
  
  real beta0;
  real<lower = 0, upper = 1> rho;
  real<lower = 0, upper = 1> alpha;
  
  // real<lower = 0, upper = 1> p_severe_k[Nk];
    
}
transformed data{
	matrix[Nk, Nk] K_L = diag_pre_multiply(sigma_k, G_t_L);
    matrix[Nt, Nt] T_L = AR1_L(Nt, rho);
    matrix[Ni, Ni] S_L = cholesky_decompose(inverse_spd(D - alpha*W));
        

}
parameters {    
    real eta[Ni, Nk, Nt];
	real<lower = 0, upper = 1> p_severe_k[Nk];
}
transformed parameters {
    real lambda_ikt[Ni, Nk, Nt];
    
    matrix[Ni, Nt] lambda_it;
    matrix[Ni, Nt] lambda_it_severe;
    
    
    
    // local variables
    {
        real theta_ikt[Ni, Nk, Nt] = rep_array(0.0,Ni, Nk, Nt); 
        
        for(t in 1:Nt)
            theta_ikt[,,t] = to_array_2d(S_L * to_matrix(eta[,,t]) * K_L');
        
        lambda_it = rep_matrix(0, Ni, Nt);
        lambda_it_severe = rep_matrix(0, Ni, Nt);
                
        for(k in 1:Nk){
            theta_ikt[,k,] = to_array_2d(to_matrix(eta[,k,]) * T_L'  + beta0);
            //lambda_ikt[,k,] = to_array_2d(exp(to_matrix(theta_ikt[,k,]) + beta0[k]));
			lambda_ikt[,k,] = to_array_2d(exp(to_matrix(theta_ikt[,k,])));
			lambda_it += to_matrix(lambda_ikt[,k,]);
            lambda_it_severe += to_matrix(lambda_ikt[,k,]) * p_severe_k[k];
        }
    }
}
model {
    to_array_1d(eta) ~ normal(0,1);
	p_severe_k ~ beta(1.1,10.9);
    
    //see https://mc-stan.org/docs/2_21/functions-reference/mixed-operations.html, when the input of to_array_1d is an array, it's converted in a row-major order, while when the input is a matrix, it's in a column-major order
    to_array_1d(Y) ~ poisson(to_array_1d((lambda_it .* N_it)'));
    to_array_1d(Y_severe) ~ poisson(to_array_1d((lambda_it_severe .* N_it)'));
    
    
    for(k in 1:Nk){
        to_array_1d(Z_mild[,k,]) ~ poisson(to_array_1d(((1-p_severe_k[k]) * to_matrix(lambda_ikt[,k,]) .* N_it .* p_test_mild)'));
        to_array_1d(Z_severe[,k,]) ~ poisson(to_array_1d((p_severe_k[k] * to_matrix(lambda_ikt[,k,]) .* N_it .* p_test_severe)'));
    }  
       
}


