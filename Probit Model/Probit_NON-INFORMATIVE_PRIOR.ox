#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl y, X, mode, H;

/*cov_matrix(const a){
  decl m, x, k, cov;
  m = meanc(a);
  x = a - m;
  k = rows(a);
  cov = x' * x / (k-1);
  return(cov);
} */

densmn(const x, const mu, const sigma2, const islog){// density of multi-normal
  decl sigma, std, logd;
  sigma = choleski(sigma2);
  std = (x - mu) * invert(sigma');
  logd = sumr(log(densn(std)));
  if(islog == 1) return(logd);
  else return(exp(logd));
}
delete_column(const x, const i){//Delete the specified column
  decl p, nx;
  p = columns(x)-1;
  if(i == 0) nx = x[][1:p];
  else{
    if(i != p) nx = x[][0:(i-1)] ~ x[][(i+1):p];
	else nx = x[][0:p-1];
  }
  return(nx);
}
probitnoinflpost(const beta, const y, const X){
  decl dbeta, n, k, pll, lF1, lF2;
  dbeta = beta';
  n = columns(dbeta);
  k = rows(dbeta);
  pll = zeros(n, 1);
  for(decl i = 0; i < n; i++){
    lF1 = log(probn(X*dbeta[][i]));
	lF2 = log(probn(-X*dbeta[][i]));
	pll[i] = sumc(y.*lF1 + (1-y).*lF2) - (2*k-1)/4 * log(dbeta[][i]' * X' * X * dbeta[][i]) + loggamma((2*k-1)/4) - (k/2)*log(M_PI);
  }
  return(pll);
}

lm(const dX, const adFunc, const avScore, const amHess){
  decl beta, n, l1, l2;
  beta = dX;
  l1 = log(probn(X*beta));
  l2 = log(probn(-X*beta));
  adFunc[0] = sumc(y.*l1 + (1-y).*l2);
  return(1);
}

hmnoinfprobit(const niter, const y, const X, const scale){
  decl p, n, beta, sigma, tildebeta, llr;
  p = columns(X);
  n = rows(X);
  beta = zeros(niter, p);
  beta[0][]	= mode';
  sigma = choleski(invert(-H));
  for(decl i = 1; i < niter; i++){
    tildebeta =  rann(1, p) * sigma' * sqrt(scale) + beta[i-1][];
	llr = probitnoinflpost(tildebeta, y, X) - probitnoinflpost(beta[i-1][], y, X);
	if(ranu(1, 1) <= exp(llr)) beta[i][] = tildebeta;
	else beta[i][] = beta[i-1][];
  }
  decl burnin = niter * 0.1;
  beta = beta[(burnin - 1):(niter - 1)][];
  return(beta);
}

main(){
  decl bank, n, p, sample, dlm;
  bank = loadmat("bank.csv");
  X = bank[][0:3];
  y = bank[][4];
  n = rows(y);
  p = columns(X);
  mode = zeros(p, 1);
  MaxBFGS(lm, &mode, &dlm, 0, 1);
  Num2Derivative(lm, mode, &H);
  sample = hmnoinfprobit(10000, y, X, 1);
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"beta1","beta2", "beta3", "beta4"});
  out.Report();
//Bayes_Factors_Calculation
  decl betak0, mk, mk0, vk, vk0, simk, simk0, usk, usk0, prob, prob0, p0; 
  mk = meanc(sample);
  vk = variance(sample);
  simk = rann(100000, p) * choleski(2 * vk)' + mk;
  usk = probitnoinflpost(simk, y, X) - densmn(simk, mk, 2*vk, 1);
  prob = log10(meanc(exp(usk)));
  
   
  p0 = p - 1; 
  prob0 = zeros(p, 1); 
  for(decl i = 0; i < p; i++){
    X = delete_column(X, i);
	mode = zeros(p0, 1);
    MaxBFGS(lm, &mode, &dlm, 0, 1);
    Num2Derivative(lm, mode, &H);
    betak0 = hmnoinfprobit(10000, y, X, 1);
	mk0 = meanc(betak0);
    vk0 = variance(betak0);
	simk0 = rann(100000, p0) * choleski(2 * vk0)' + mk0;
	usk0 = probitnoinflpost(simk0, y, X) - densmn(simk0, mk0, 2*vk0, 1);
	prob0[i] = log10(meanc(exp(usk0)));
	
	X = bank[][0:3];
  }	
  println("%r",{"beta1", "beta2", "beta3", "beta4"},
	        "%c",{"Bayes Factors"},
	        (prob - prob0));
}