#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl X, y, mode, H;

loglinll(const dbeta, const y, const X){
  decl lF, pll, beta, n;
  beta = dbeta';
  n = columns(beta);
  pll = zeros(n, 1);
  for(decl i = 0; i < n; i++){
	lF = exp(X*beta[][i]);
    pll[i] = sumc(log(denspoisson(y, lF)));
  }
  return(pll);
}
lm(const dX, const adFunc, const avScore, const amHess){
  decl lF, beta;
  beta = dX;
  lF = exp(X*beta);
  adFunc[0] = sumc(log(denspoisson(y, lF)));
  return(1);
}
hmflatloglin(const niter, const y, const X, const scale){
  decl p, n, sigma, beta, tildebeta, llr;
  p = columns(X);
  n = rows(X);
  beta = zeros(niter, p);
  beta[0][] = mode';
  sigma = choleski(invert(-H));
  for(decl i = 1; i < niter; i++){
    tildebeta =  rann(1, p) * sigma' * sqrt(scale) + beta[i-1][];
	llr = loglinll(tildebeta, y, X) - loglinll(beta[i-1][], y, X);
	if(ranu(1, 1) <= exp(llr)) beta[i][] = tildebeta;
	else beta[i][] = beta[i-1][];
  }
  decl burnin = niter * 0.1;
  beta = beta[burnin:(niter - 1)][];
  return(beta);
}

main(){
  decl airquality, n, p, sample, dlm;
  airquality = loadmat("airquality_after_treatment.csv");
  y = airquality[][0];
  X = airquality[][1:];
  n = rows(y);
  p = columns(X);
  mode = 0.1*ones(p, 1);
  MaxBFGS(lm, &mode, &dlm, 0, 1);
  Num2Derivative(lm, mode, &H);
  sample = hmflatloglin(10000, y, X, 0.5);
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"beta1","beta2", "beta3", "beta4",
                   "beta5","beta6", "beta7", "beta8",
				   "beta9","beta10", "beta11", "beta12",
				   "beta13","beta14", "beta15", "beta16"});
  out.Report();
}