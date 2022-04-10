#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl X, y, mode, H;

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

loglinnoinflpost(const dbeta, const y, const X){
  decl lF, pll, beta, n, k;
  beta = dbeta';
  n = columns(beta);
  k	= rows(beta);
  pll = zeros(n, 1);
  for(decl i = 0; i < n; i++){
	lF = exp(X*beta[][i]);
    pll[i] = sumc(log(denspoisson(y, lF))) - k/2 * log(beta[][i]' * X' * X * beta[][i]/2) + loggamma(k/2) - (k/2)*log(M_2PI);
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
hmnoinfloglin(const niter, const y, const X, const scale){
  decl p, n, sigma, beta, tildebeta, llr;
  p = columns(X);
  n = rows(X);
  beta = zeros(niter, p);
  beta[0][] = mode';
  sigma = choleski(invert(-H));
  for(decl i = 1; i < niter; i++){
    tildebeta =  rann(1, p) * sigma' * sqrt(scale) + beta[i-1][];
	llr = loglinnoinflpost(tildebeta, y, X) - loglinnoinflpost(beta[i-1][], y, X);
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
  sample = hmnoinfloglin(10000, y, X, 0.5);
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"beta1","beta2", "beta3", "beta4",
                   "beta5","beta6", "beta7", "beta8",
				   "beta9","beta10", "beta11", "beta12",
				   "beta13","beta14", "beta15", "beta16"});
  out.Report();
//Bayes_Factors_Calculation
  decl betak0, mk, mk0, vk, vk0, simk, simk0, usk, usk0, prob, prob0, p0; 
  mk = meanc(sample);
  vk = variance(sample);
  simk = rann(100000, p) * choleski(2 * vk)' + mk;
  usk = loglinnoinflpost(simk, y, X) - densmn(simk, mk, 2*vk, 1);
  prob = log10(meanc(exp(usk)));
   
  p0 = p - 1; 
  prob0 = zeros(p, 1); 
  for(decl i = 0; i < p; i++){
    X = delete_column(X, i);
	mode = zeros(p0, 1);
    MaxBFGS(lm, &mode, &dlm, 0, 1);
    Num2Derivative(lm, mode, &H);
    betak0 = hmnoinfloglin(10000, y, X, 1);
	mk0 = meanc(betak0);
    vk0 = variance(betak0);
	simk0 = rann(100000, p0) * choleski(2 * vk0)' + mk0;
	usk0 = loglinnoinflpost(simk0, y, X) - densmn(simk0, mk0, 2*vk0, 1);
	prob0[i] = log10(meanc(exp(usk0)));
	
	X = airquality[][1:];
  }	
  println("%r",{"beta1","beta2", "beta3", "beta4",
                "beta5","beta6", "beta7", "beta8",
				"beta9","beta10", "beta11", "beta12",
				"beta13","beta14", "beta15", "beta16"},
	        "%c",{"Bayes Factors"},
	        (prob - prob0));
}