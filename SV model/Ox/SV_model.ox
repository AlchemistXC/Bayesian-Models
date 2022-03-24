#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl y, n;
decl mu0, sm0, sm02, nu0, delta0, a0, b0;
decl phi, sig2, sig, mu, h;
decl t;
decl modeh, mlh, Hh, vh2, vh, mh, sh;

progress_bar(const i, const iter, const burnin);
ranmu(const dphi, const dsig2, const dh);
ransigma2(const dphi, const dmu, const dh);
rantrancenorm(const r, const c, const mu, const sigma, const a, const b);
ranphiprop(const dsig2, const dmu, const dh);
flogphi(const dphi, const dsig2, const dmu, const dh);
ranhprop(const dphi, const dmu, const dsig2, const dh);
flogh(const dht, const dphi, const dmu, const dsig2, const dh);

main(){
  decl iter, burnin, irs, sample, sampleh;
  decl ach, acphi;
  decl rho, u, p, q, prop;

  // Data
  y = loadmat("USCPI.csv");
  y = y - meanc(y);
  n = rows(y);

// Set initial value
  phi = 0.5;	//phi
  sig = 0.5;	//sigma
  sig2 = sig^2;	//sigma^2
  mu  = 1;	//mu  
  h = zeros(n, 1); //h
// Set prior
  a0 = 0.5;  //(phi+1)/2 ~ Beta(ap0, bp0)
  b0 = 0.5;
  nu0 = 0.01;  //sig^2 ~ InvGamma(n0/2, S0/2)
  delta0 = 0.01;
  mu0 = 0;  //mu ~ N(mm0, sm0^2)
  sm0 = 1;
  sm02  = sm0^2;
// Set MCMC
  iter = 5000;
  burnin = 0.1*iter;
  irs = 7;
  ranseed(irs);
  sample = zeros(iter, 3);
  sampleh = zeros(iter, n);

// Set MH
  ach = 0;	// all acceptance count of h
  acphi = 0;  //acceptance count of phi
  for(decl i = -burnin; i < iter; i++){
	progress_bar(i, iter, burnin);
	// Simulate mu
	mu = ranmu(phi, sig2, h);
	// Simulate sigma2
	sig2 = ransigma2(phi, mu, h);
	// Simulate phi (MH)
	prop = ranphiprop(sig2, mu, h);
	p =   flogphi(prop[0], sig2, mu, h)
	    + 0.5 * (prop[0] - prop[1])^2 / prop[2];
	q =   flogphi(phi, sig2, mu, h)
	    + 0.5 * (phi - prop[1])^2 / prop[2];
	rho = exp(p - q);
	  u = ranu(1, 1);
	  if(u <= rho){
	    phi = prop[0];
		acphi++;
	  }
	//Simulate h (MH)
    for(t = 0; t < n; t++){
	  prop = ranhprop(phi, mu, sig2, h);
	  p =  flogh(prop[0], phi, mu, sig2, h)
	     + 0.5 * (prop[0] - prop[1])^2 / prop[2];
	  q =  flogh(h[t], phi, mu, sig2, h)
	     + 0.5 * (h[t] - prop[1])^2 / prop[2];
	  rho = exp(p - q);
	  u = ranu(1, 1);
	  if(u <= rho){
	    h[t] = prop[0];
		ach++;
	  }
	}
	//Record data when i>=0
	if(i >= 0){
	  sample[i][] = mu ~ phi ~ sqrt(sig2);
      sampleh[i][] = h';
	}
  }	
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"mu","phi","sigma"});
  out.Report();
  decl meanh = meanc(sampleh);
  Draw(0, meanh);
  SaveDrawWindow("h.eps");
  println("Acceptance rate of phi (MH): ", acphi/(burnin+iter)*100, "%");
  println("Average acceptance rate of each element in vector h (MH): ", ach/((burnin+iter)*n)*100,"%");
}// End of main
//-----------------------------------------------


progress_bar(const i, const iter, const burnin){
  decl form1, form2, fnum;
  form1 = round(<1:9>*(iter-1)/10);
  form2 = <1:9>*10;
  fnum = 0;
  if(i == -burnin){println("Burn-in period is in progress.");}
  else{
	if(i == 0){println("Burn-in period is complete.");
	           println("Samples are being collected.");
			   print("Accomplishment rate: 0% -> ");}
	else{
	  do{fnum++;}while(i > form1[fnum] && fnum < 8);
	  if(i == form1[fnum]){print(form2[fnum],"% -> ");}
	  else{
	    if(i == round((iter-1))){println("100%");
	                             println("All phases complete.");}
	  }
	}		   
  }
}

ranmu(const dphi, const dsig2, const dh){
  decl c1, c2, c3, mu1, sm12, sm1, dmu;
  c1 = (n-1)*(1-dphi)^2 + 1 - dphi^2;
  sm12 = sm02 * dsig2 /(sm02 * c1 + dsig2);
  sm1 = sqrt(sm12);
  c2 = sumc(dh[1:(n-1)] - dphi*dh[0:(n-2)]);
  c3 = ((1 - dphi^2)*dh[0] + (1 - dphi)*c2) / dsig2;
  mu1 = sm12 * (c3 + mu0/sm02);
  dmu = sm1 * rann(1, 1) + mu1;
  return(dmu);
}

ransigma2(const dphi, const dmu, const dh){
  decl nu1, delta1, c1, c2, dsig2;
  nu1 = n + nu0;
  c1 = dh - dmu;
  c2 = sumc((c1[1:(n-1)] - dphi*c1[0:(n-2)]).^2);
  delta1 = delta0 + (dh[0] - dmu)^2 * (1-dphi^2) + c2;
  dsig2 = 1 / rangamma(1, 1, nu1/2, delta1/2);
  return(dsig2);
}

rantrancenorm(const r, const c, const mu, const sigma, const a, const b){// tranced normal distribution
  decl x, low, high, za, zb, u;
  za = (a-mu)/sigma;
  zb = (b-mu)/sigma;
  low = probn(za);
  high = probn(zb);
  u = ranu(r, c) * (high - low) + low;
  x = quann(u);
  x = sigma * x + mu;
  return(x);
}

ranphiprop(const dsig2, const dmu, const dh){
  decl muphi, sig2phi, sigphi, c1, c2, c3, dphi;
  c1 = dh - dmu;
  c2 = sumc(c1[0:(n-2)].^2);
  c3 = sumc(c1[1:(n-1)] .* c1[0:(n-2)]);
  sig2phi = dsig2 / c2;
  sigphi = sqrt(sig2phi);
  muphi = c3 / c2;
  dphi = rantrancenorm(1, 1, muphi, sigphi, -1, 1);
  return(dphi ~ muphi ~ sig2phi);
}

flogphi(const dphi, const dsig2, const dmu, const dh){
  decl c1, c2, c3, c4, logf;
  c1 = 1 - dphi^2;
  c2 = dh - dmu;
  c3 = sumc((c2[1:(n-1)] - dphi*c2[0:(n-2)]).^2); 
  logf =   log(densbeta((dphi+1)/2,a0,b0))
         + 0.5 * log(c1)
		 - 0.5 * c3 / dsig2
		 - 0.5 * (c2[0].^2) * c1 / dsig2;
  return(logf);		 
}

ranhprop(const dphi, const dmu, const dsig2, const dh){
  decl sig2t, ht, mut, c, dht;
  if(t == 0){
    sig2t = dsig2;
	ht = dphi * ((1 - dphi) * dmu + dh[1]); 
  }
  else{
    if(t != n-1){
	  c = 1 + dphi^2;
	  sig2t = dsig2 / c;
	  ht = dphi * (dh[t+1] + dh[t-1])/ c;  
	}
	else{
	  sig2t = dsig2;
	  ht = dphi * dh[n-2] + (1 - dphi) * dmu;
	}
  }
  mut = ht ;//+ 0.5 * sig2t * ((y[t]^2) * exp(-ht) - 1);
  dht = sqrt(sig2t) * rann(1, 1) + mut;
  return(dht ~ mut ~ sig2t);
}

flogh(const dht, const dphi, const dmu, const dsig2, const dh){
  decl sig2t, ht, c, f;
  f = - 0.5 * (y[t]^2 * exp(-dht) + dht);
  if(t == 0){
    sig2t = dsig2;
	ht = dphi * ((1 - dphi) * dmu + dh[1]); 
  }
  else{
    if(t != n-1){
	  c = 1 + dphi^2;
	  sig2t = dsig2 / c;
	  ht = dphi * (dh[t+1] + dh[t-1])/ c;  
	}
	else{
	  sig2t = dsig2;
	  ht = dphi * dh[n-2] + (1 - dphi) * dmu;
	}
  }
  f = f - 0.5 * (dht - ht)^2 / sig2t;
  return(f);
}