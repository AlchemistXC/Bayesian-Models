#include<oxstd.h>
#include <oxprob.h>
#include <ReportMCMC.ox>

main(){
  decl data, y, x, n, k;
  decl beta0, B0, inv_B0, alpha0, delta0;
  decl beta1, B1, alpha1, delta1;
  decl beta, sig, phi;
  decl iter, burnin, rs, sample;
// Loading data  
  data = loadmat("iris.csv");
  n = rows(data);
  y = data[][0];
  x = ones(n, 1) ~ data[][1:3];
  k = columns(x);
// prior
  beta0 = zeros(k, 1);
  B0 = diag(ones(k, 1));
  inv_B0 = invert(B0);
  alpha0 = 0.001;
  alpha1 = alpha0 + n;//posterior
  delta0 = 0.001;
// Initial value
  beta = ones(k, 1);
  sig = 1;
// MCMC set
  iter = 5000;
  burnin = 0.1 * iter;
  rs = 7;
  ranseed(rs);
  sample = zeros(iter, k+1);
// MCMC
  for(decl i = -burnin; i < iter; i++){
	//Simulate sigma2;
	delta1 = (y - x * beta)' * (y - x * beta) + delta0;
	phi = rangamma(1, 1, alpha1/2, delta1/2);
	sig = sqrt(1 / phi);
	//Simulate beta;
	B1 = invert(phi * x' * x + inv_B0);
	beta1 = B1 * (phi * x' * y + inv_B0 * beta0);
	beta = choleski(B1) * rann(k, 1) + beta1;
	// Sample
	if(i >= 0){
	  sample[i][] = sig ~ beta';
	}
  }
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"sigma","Beta0", "Beta1", "Beta2", "Beta3"});
  out.Report();
}