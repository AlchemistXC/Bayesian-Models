#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl X, y, mode, H;
logitll(const dbeta, const y, const X){
  decl lF1, lF2, p, beta;
  beta = dbeta';
  lF1 = - log(1 + exp(-X*beta));
  lF2 = - log(1 + exp(X*beta));
  p = sumc(y .* lF1 + (1 - y) .* lF2);
  return(p);
}
lm(const dX, const adFunc, const avScore, const amHess){
  decl lF1, lF2, beta;
  beta = dX;
  lF1 = - log(1 + exp(-X*beta));
  lF2 = - log(1 + exp(X*beta));
  adFunc[0] = sumc(y .* lF1 + (1 - y) .* lF2);
  return(1);
}
hmflatlogit(const niter, const y, const X, const scale){
  decl p, n, sigma, beta, tildebeta, llr;
  p = columns(X);
  n = rows(X);
  beta = zeros(niter, p);
  beta[0][] = mode';
  sigma = choleski(invert(-H));
  for(decl i = 1; i < niter; i++){
    tildebeta =  rann(1, p) * sigma' * sqrt(scale) + beta[i-1][];
	llr = logitll(tildebeta, y, X) - logitll(beta[i-1][], y, X);
	if(ranu(1, 1) <= exp(llr)) beta[i][] = tildebeta;
	else beta[i][] = beta[i-1][];
  }
  decl burnin = niter * 0.1;
  beta = beta[burnin:(niter - 1)][];
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
  sample = hmflatlogit(10000, y, X, 1);
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"beta1","beta2", "beta3", "beta4"});
  out.Report();	 
}