#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

dnorm(const x, const mean, const sd, const islog){
  decl d;
  d = sqrt(0.5*M_1_PI) ./ sd .* exp(-0.5 * ((x - mean).^2) ./ (sd.^2)); 
  if(islog == 0) return(d);
  else return(log(d));
}

sum(x){
  if(rows(x) == 1) return(sumr(x));
  else{
    if(columns(x) == 1) return(sumc(x));
	else return(sumr(sumc(x)));
  }
}

plotmix(const n, const p, const dmu1, const dsigma1, const dmu2, const dsigma2, const choose){
  decl n1, pbar, mu1, mu2, dif, length, mo1, mo2, ca1, ca2, like, sampl;
  n1 = round(n*p);
  sampl = (rann(n1, 1) * dsigma1 + dmu1) | (rann(n-n1, 1) * dsigma2 + dmu2);
  pbar = 1 - p;
  dif = maxc(sampl) - minc(sampl);
  length = round(dif / 0.1);
  mu1 = zeros(length, 1);
  mu1[0][] = minc(sampl);
  for(decl i = 1; i < length; i++){
    mu1[i][] = mu1[i-1][] + 0.1;
  }
  mu2 = mu1;
  mo1 = mu1 * ones(1, length);
  mo2 = ones(length, 1) * mu2';
  ca1 = -0.5 * mo1 * mo1;
  ca2 = -0.5 * mo2 * mo2;
  like = 0*mo1;
  for(decl t = 0; t < n; t++){
    like += log(p * exp(ca1 + sampl[t]*mo1)+pbar*exp(ca2+sampl[t]*mo2));
  }
  like += 0.1*(ca1 + ca2);
  if(choose == 1) return(like);
  else return(sampl);
}

lpost(const x, const mu, const p, const delta, const lambda){
  decl l;
  l = sum(log(p * dnorm(x, mu[0], 1, 0) + (1-p) * dnorm(x,mu[1], 1, 0))) +
      sum(dnorm(mu, delta, 1/sqrt(lambda), 1));
  return(l);	  
}

hmmeantemp(const dat, const niter, const var, const alpha){
  decl mu, muprop, bound;
  mu = zeros(niter, 2);
  mu[0][] = <1, 3>;
  for(decl i = 1; i < niter; i++){
    muprop = rann(1, 2)*sqrt(var)' + mu[i-1][];
	bound = lpost(dat, muprop, 0.7, 0, 1) - lpost(dat, mu[i-1][], 0.7, 0, 1);
	if(ranu(1, 1) <= exp(alpha * bound)) mu[i][] = muprop;
	else mu[i][] = mu[i-1][];
  }
  return(mu);
}

main(){
  decl dat, simu;
  dat = plotmix(10^4, 0.7, 2.5, 1, 0, 1, 0);
  simu = hmmeantemp(dat,10^4, 1, 1);
  //simu = hmmeantemp(dat,10^4, 1, 0.1);
  //simu = hmmeantemp(dat,10^4, 1, 0.01);
  decl out = new ReportMCMC(simu);
  out.SetVarNames({"mu1","mu2"});
  out.Report();
}