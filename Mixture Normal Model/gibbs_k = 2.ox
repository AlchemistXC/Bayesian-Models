#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>



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

omega(const z, const x, const p){
  decl n, n1, n2, xbar1, xbar2, z1, z2, ss1, ss2;
  n = rows(z);
  n1 = 0;
  n2 = 0;
  z1 = zeros(n, 1);
  z2 = zeros(n, 1);
  for(decl i = 0; i < n; i++){
    if(z[i] == 1){
	   n1++;
	   z1[i] = 1;
	}
	if(z[i] == 2){
	  n2++;
	  z2[i] = 1;
	}
  }
  if(n1 == 0) xbar1 = 0;
  else xbar1 = (z1' * x)/n1;
  if(n2 == 0) xbar2 = 0;
  else xbar2 = (z2' * x)/n2;
  ss1 = z1' * (x - xbar1).^2;
  ss2 = z2' * (x - xbar2).^2;
  return (sqrt((n1 + .25) * (n2 + .25)) *
          p^n1 * (1-p)^n2 *
		  exp(-((n1 + .25)*ss1 + (n2 + .25)*ss2)/2) *
		  exp(-(n1*xbar1^2 + n2*xbar2)/8)
		  );  
}

gibbsmean(const p, const datha, const niter){
  decl n, z, ssiz, nxj, mug, prob, nj, zj;
  n = rows(datha);
  z = zeros(n, 1);
  ssiz = zeros(1, 2);
  nxj = zeros(1, 2);
  prob = zeros(1, 2);
  mug = meanc(datha) * ones(niter+1, 2);
  for(decl i = 1; i < niter+1; i++){
    nj = zeros(2, 1);
	zj = zeros(n, 2);
	for(decl t = 0; t < n; t++){
	  prob[0] = p * exp(-0.5*(datha[t]-mug[i-1][0])^2);
	  prob[1] = (1-p) * exp(-0.5*(datha[t]-mug[i-1][1])^2);
	  z[t] = ranbinomial(1, 1, 1, prob[1]/sumr(prob)) + 1;
	  if(z[t] == 1){
		nj[0]++;
		zj[t][0] = 1;
	  }
	  if(z[t] == 2){
		nj[1]++;
		zj[t][1] = 1;
	  }	
	}
	for(decl j = 0; j < 2; j++){
	  ssiz[j] = 1 + nj[j];
	  nxj[j] = zj[][j]' * datha;
	}
	mug[i][] = rann(1, 2) .* sqrt(1/ssiz) + (nxj ./ ssiz);
  }
  return(mug);
}
main(){
  decl like, sampl, z, fit, dat;
  sampl = plotmix(4, 0.7, 2.5, 1, 0, 1, 0);
  z = ranbinomial(4, 1, 1, 0.5) + 1;
  println("Omega:", omega(z, sampl, 0.8));
  dat =	plotmix(500, 0.7, 2.5, 1, 0, 1, 0);
  fit = gibbsmean(0.7, dat, 10^4);
  decl out = new ReportMCMC(fit);
  out.SetVarNames({"mu1","mu2"});
  out.Report();
}