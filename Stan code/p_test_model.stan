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
  int<lower=1> Nt;   // number of time point
  
  int Y[Ni, Nt];      // number of all cases
  int Z[Ni, Nt];  //number of tested cases
}
parameters {    
    
    matrix[Ni, Nt] eta;
    
    real<lower = 0, upper = 1> rho; //temporal autocorrelation XXfor each prefecture
    
    real beta0; //overall probability
	real beta_pref[Ni]; //prefecture random effects
}
transformed parameters {
    matrix[Ni, Nt] theta_it;
	
	//local variables
	{
		matrix[Nt, Nt] T_L = AR1_L(Nt, rho);	
		matrix[Ni, Nt] eta_transform;
		
		// for fast computation, we generated correlated random numbers with uncorrelated random numbers and the cholesky decomposition of the covariance matrix, see https://math.stackexchange.com/questions/163470/generating-correlated-random-numbers-why-does-cholesky-decomposition-work
		eta_transform = eta * T_L';
		
		for(i in 1:Ni){
			
		theta_it[i, ] = inv_logit(beta0 + beta_pref[i] + eta_transform[i,]);
		}
	}
}
model {
    rho ~ beta(2,2);
    to_vector(eta) ~ normal(0,1);
    
    beta0 ~ student_t(3,0,1);
	beta_pref ~ normal(0,1);
    
	to_array_1d(Z) ~ binomial(to_array_1d(Y), to_array_1d(theta_it'));
}
