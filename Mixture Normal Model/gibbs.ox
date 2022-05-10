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
decl mug, sigg, prog, lopost;
gibbsnorm(const datha, const k, const niter){
  decl n, mu, sig, z, zj, nxj, ssiz, ssum, lik, prob, i, t, j, nu, s, epsilon, sig_mu;
  n = rows(datha);
  mu = meanc(datha);
  sig = varc(datha);
  z = zeros(n, 1);
  zj = zeros(n, k);
  nxj = zeros(1, k);
  ssiz = zeros(1, k);
  ssum = zeros(1, k);
  mug = zeros(niter, k);
  sigg = zeros(niter, k);
  prog = zeros(niter, k);
  lopost = zeros(niter, 1);
  lik = zeros(n, k);
  nu = zeros(k, 1);
  s = zeros(k, 1);
  epsilon =	zeros(k, 1);
  sig_mu = zeros(k, 1);
  prog[0][] = ones(1, k) / k;
  mug[0][] = mu * ones(1, k);
  sigg[0][]	= sig * ones(1, k);
  // current log-likelihood
  for(decl j = 0; j < k; j++){
    //lik[][j] = prog[0][j] * sqrt(sigg[0][j]) * densn((datha-mug[0][j])/sqrt(sigg[0][j]));
	lik[][j] = prog[0][j] * dnorm(datha, mug[0][j], sqrt(sigg[0][j]), 0);
	
  }

  lopost[0] = sum(log(sumr(lik))) +
	          sum(dnorm(mug[0][], mu, sqrt(sigg[0][]), 1)) -
			  11 * sum(log(sigg[0][])) -
			  sum(sig./sigg[0][]) +
			  0.5 * sum(log(prog[0][]));
			  
  for(i = 0; i < (niter-1); i++){
    for(t = 0; t < n; t++){
	//missing data completion
      prob = prog[i][] .* dnorm(datha[t], mug[i][], sqrt(sigg[i][]), 0);
	  if(sum(prob) == 0) prob = ones(1, k) / k;
	  zj[t][] =	 ranmultinomial(1, prob'/sum(prob));
	  z[t] = zj[t][] * <1, 2, 3>'; 
	}
	for(j = 0; j < k; j++){
	  ssiz[j] = sum(zj[][j]);
	  nxj[j] = zj[][j]' * datha;
	}
	epsilon = (mu + nxj) ./ (ssiz + 1);
	sig_mu = sigg[i][] ./ (ssiz + 1);
	mug[i+1][] = rann(1, k) .* sqrt(sig_mu) + epsilon;
	for(j = 0; j < k; j++){
	  ssum[j] = zj[][j]' * (datha - nxj[j]/ssiz[j]).^2;
	  nu[j] = 20 + ssiz[j];
	  s[j] = 2 * sig + ssum[j] + ssiz[j]/(ssiz[j] + 1) * (epsilon[j] - zj[][j]'*datha / ssiz[j])^2;
	  sigg[i+1][j] = 1 / rangamma(1, 1, nu[j]/2, s[j]/2);
	}

	prog[i+1][0:(k-2)] = randirichlet(1, ssiz + 0.5); //different from R
	prog[i+1][k-1] = 1 - sum(prog[i+1][0:(k-2)]);
    //current log-likelihood
    for(j = 0; j < k; j++){
	  lik[][j] = prog[i+1][j] * dnorm(datha, mug[i+1][j], sqrt(sigg[i+1][j]), 0);
    }
	lopost[i+1] = sum(log(sumr(lik))) +
	              sum(dnorm(mug[i+1][], mu, sqrt(sigg[i+1][]), 1)) -
			      11 * sum(log(sigg[i+1][])) -
			      sum(sig./sigg[i+1][]) +
			      0.5 * sum(log(prog[i+1][]));
  }
  return 1;
}

main(){
  decl datha, n, mu, sig, fit;
  datha = loadmat("datha.csv");
  n = rows(datha);
  mu = meanc(datha);
  sig = varc(datha)*n/(n-1);
  gibbsnorm(datha, 3, 1000);
  
// Label
  decl indimap, m_mu, m_p, m_sig;
  indimap = maxcindex(lopost);
  m_mu = mug[indimap][];
  m_p = prog[indimap][];
  m_sig = sigg[indimap][];
  println("mu:", m_mu);
  println("sig:", m_sig);
  println("p:", m_p);
  decl lili, alloc, ormu, orsig, orp, perma, t, j, i, entropies, best;
  lili = alloc = zeros(n, 3);
  for(t = 0; t < n; t++){
	lili[t][] = m_p .* dnorm(datha[t], m_mu, sqrt(m_sig), 0); 
	lili[t][] = lili[t][] / sum(lili[t][]);
  }
  ormu = orsig = orp = zeros(1000, 3);
  perma = <1, 2, 3;
           1, 3, 2;
		   3, 1, 2;
		   3, 2, 1;
		   2, 1, 3;
		   2, 3, 1> - 1;// k = 3;
  for(t = 0; t < 1000; t++){
    entropies = zeros(6, 1);
	for(j = 0; j < n; j++){
	  alloc[j][] = prog[t][] .* dnorm(datha[j], mug[t][], sqrt(sigg[t][]), 0);
	  alloc[j][] = alloc[j][] / sum(alloc[j][]);
	  for(i = 0; i < 6; i++){
	    entropies[i] += sum(lili[j][] .* log(alloc[j][perma[i][]]));
	  }
	}
	best = maxcindex(entropies);
	ormu[t][] = mug[t][perma[best][]];
	orsig[t][] = sigg[t][perma[best][]];
	orp[t][] = prog[t][perma[best][]];
  }
  decl out = new ReportMCMC(ormu);
  out.SetVarNames({"mu1", "mu2", "mu3"});
  out.Report();
}