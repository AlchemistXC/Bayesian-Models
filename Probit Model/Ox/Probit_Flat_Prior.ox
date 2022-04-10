#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl y, X, mode, cov;
probitll(const beta, const y, const X){
  decl dbeta, n, pll, lF1, lF2;
  dbeta = beta';
  n = columns(dbeta);
  pll = zeros(n, 1);
  for(decl i = 0; i < n; i++){
    lF1 = log(probn(X*dbeta[][i]));
	lF2 = log(probn(-X*dbeta[][i]));
	pll[i] = sumc(y.*lF1 + (1-y).*lF2);
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
hmflatprobit(const niter, const y, const X, const scale){
  decl p, n, beta, sigma, tildebeta, llr;
  p = columns(X);
  n = rows(X);
  beta = zeros(niter, p);
  beta[0][]	= mode';
  sigma = choleski(cov);
  for(decl i = 1; i < niter; i++){
    tildebeta =  rann(1, p) * sigma' * sqrt(scale) + beta[i-1][];
	llr = probitll(tildebeta, y, X) - probitll(beta[i-1][], y, X);
	if(ranu(1, 1) <= exp(llr)) beta[i][] = tildebeta;
	else beta[i][] = beta[i-1][];
  }
  return(beta);
}

main(){
  decl bank, n, p, dlm, H, sample;
  bank = loadmat("bank.csv");
  X = bank[][0:3];
  y = bank[][4];
  n = rows(y);
  p = columns(X);
  mode = zeros(p, 1);
  MaxBFGS(lm, &mode, &dlm, 0, 1);
  Num2Derivative(lm, mode, &H);
  cov = -invert(H);
  sample = hmflatprobit(10000, y, X, 1)[1000:9999][];
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"beta1","beta2", "beta3", "beta4"});
  out.Report();
}