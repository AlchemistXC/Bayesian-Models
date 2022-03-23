#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>
#include <ReportMCMC.ox>

decl y, n;
decl mu0, sm0, sm02, nu0, delta0, a0, b0;
decl t, isp;
decl phi, sig2, sig, mu, h;
decl modephi, modeh, mlphi, mlh, Hphi, Hh, vphi, vh, mphi, mh, sphi, sh;
  decl logc, ma, mx;
progress_bar(const i, const iter, const burnin);
ranmu(const dphi, const dsig2, const dh);
ransigma2(const dphi, const dmu, const dh);
flogphi(const dX, const adFunc, const avScore, const amHess);
logphi(const dphi, const dsig2, const dmu, const dh);
flogh(const dX, const adFunc, const avScore, const amHess);
logph(const dht, const dh, const dphi, const dsig2, const dmu, const dt);
rantrancenorm(const r, const c, const mu, const sigma, const a, const b);

main(){
  decl iter, burnin, irs, sample, sampleh;
  decl acphi, ach;
  decl logrhophi, logrhoh, rhophi, rhoh, logp, logpphi, ratio, prop, u;
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
  a0 = 20;  //(phi+1)/2 ~ Beta(ap0, bp0)
  b0 = 1.5;
  nu0 = 5;  //sig^2 ~ InvGamma(n0/2, S0/2)
  delta0 = 0.05;
  mu0 = 1.5;  //mu ~ N(mm0, sm0^2)
  sm0 = 0.5;
  sm02  = sm0^2;
// Set MCMC
  iter = 5000;
  burnin = 0.1*iter;
  irs = 7;
  ranseed(irs);
  sample = zeros(iter, 3);
  sampleh = zeros(iter, n);

// Set MH
  acphi = 0; // acceptance count of phi
  ach = 0;	// all acceptance count of h


// MCMC
  for(decl i = -burnin; i < iter; i++){
	progress_bar(i, iter, burnin);
	// Simulate h　(AR-MH)
	for(t = 0; t < n; t++){
	  isp = 0; //compute mode
	  MaxBFGS(flogh, &modeh, &mlh, 0, 1);
	  Num1Derivative(flogh, modeh, &sh);
	  Num2Derivative(flogh, modeh, &Hh);
	  vh = sqrt(-1/Hh);
	  mh = modeh + vh*vh * sh;
	  isp = 1; //compute c:max(fx/cgx)=1
      MaxBFGS(flogh, &mx, &ma, 0, 1);
	  logc = ma - 1;
	  do{
	    prop = vh * rann(1, 1) + mh; 
	    logp = logph(prop, h, phi, sig2, mu, t);
		u = ranu(1, 1);
	  }while(u > exp(logp));
	  ratio = logph(h[t], h, phi, sig2, mu, t);
	  if(ratio < 0){logrhoh = 0;}
      else{
		if(logp < 0){logrhoh = -ratio;}
		else{logrhoh = logp - ratio;}
	  }
	  rhoh = exp(logrhoh);
	  u = ranu(1, 1);
	  if(u <= rhoh){
		h[t] = prop;
		ach++;
	  } 
	}
	// Simulate sigma2
	sig2 = ransigma2(phi, mu, h);
	// Simulate phi	(AR-MH)
	isp = 0; //compute mode
	modephi = 0; // initial value
	MaxBFGS(flogphi, &modephi, &mlphi, 0, 1);
	Num1Derivative(flogphi, modephi, &sphi);
	Num2Derivative(flogphi, modephi, &Hphi);
	vphi = sqrt(-1/Hphi);
	mphi = modephi + vphi*vphi*sphi;
	isp = 1; //compute c:max(fx/cgx)=1
	mx = 0; // initial value
	MaxBFGS(flogphi, &mx, &ma, 0, 1);
	logc = ma - 1;
	do{
	  prop = rantrancenorm(1, 1, mphi, vphi, -1, 1);
	  logpphi = logphi(prop, sig2, mu, h);
	  u =ranu(1, 1);
	}while(u > logpphi);
	ratio = logphi(phi, sig2, mu, h);
	if(ratio < 0){logrhophi = 0;}
	else{
	  if(logpphi < 0){logrhophi = -ratio;}
	  else{logrhophi = logpphi - ratio;}
	}
	rhophi = exp(logrhophi);
	u = ranu(1, 1);
	if(u <= rhophi){
	phi = prop;
	acphi++;
	}
	// Simulate mu
	mu = ranmu(phi, sig2, h);
	//Record data when i>=0
	if(i >= 0){
	  sample[i][] = mu ~ phi ~ sqrt(sig2);
      sampleh[i][] = h';
	}
  }
  decl out = new ReportMCMC(sample);
  out.SetVarNames({"mu","phi","sigma"});
  out.Report();
  println("Acceptance rate of phi: ", acphi/(burnin+iter)*100, "%");
  println("Average acceptance rate of each element in vector h: ", ach/((burnin+iter)*n)*100,"%");
  decl meanh = meanc(sampleh);
  Draw(0, meanh);
  SaveDrawWindow("h.eps");
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

flogphi(const dX, const adFunc, const avScore, const amHess){
  decl dphi, c1, c2, c3, logf;
  dphi = dX;
  c1 = 1 - dphi^2;
  c2 = h - mu;
  c3 = sumc((c2[1:(n-1)] - dphi*c2[0:(n-2)]).^2); 
  logf =   log(0.5*densbeta((dphi+1)/2,a0,b0))
              + 0.5 * log(c1)
			  - 0.5 * c3 / sig2
			  - 0.5 * (c2[0].^2) * c1 / sig2;
  if(isp == 0){adFunc[0] = logf;}
  else{adFunc[0] = logf - log(densn((dphi-mphi)/vphi));}
  return 1;			   
}

logphi(const dphi, const dsig2, const dmu, const dh){
  decl f, c1, c2, c3;
  c1 = 1 - dphi^2;
  c2 = dh - dmu;
  c3 = sumc((c2[1:(n-1)] - dphi*c2[0:(n-2)]).^2);
  f = log(0.5*densbeta((dphi+1)/2,a0,b0))
	  + 0.5 * log(c1)
	  - 0.5 * c3 / dsig2
	  - 0.5 * (c2[0].^2) * c1 / dsig2
	  - log(densn((dphi - mphi)/vphi))
	  - logc;
  return(f);		  
}

flogh(const dX, const adFunc, const avScore, const amHess){
  decl x, c1, c2, logf;
  x = dX; 		   
  if(t == 0){
    c1 = x | h[1] - mu;
	c2 = ((c1[1] - phi * c1[0]).^2)/sig2;	
	logf = - 0.5 * x
	            - 0.5 * (y[t].^2) * exp(-x)
				- 0.5 * (1-phi^2) / sig2 * (c1[0].^2)
				- 0.5 * c2; 
  }
  else{
    if(t != n-1){
	  c1 = h[t-1] | x | h[t+1] - mu;
	  c2 = ((c1[1:2] - phi * c1[0:1]).^2)/sig2;	  
	  logf = - 0.5 * x
	              - 0.5 * (y[t].^2) * exp(-x)
				  - 0.5 * sumc(c2);
	}
	else{
	  c1 = h[t-1] | x - mu;
	  c2 = ((c1[0] - phi * c1[1]).^2)/sig2;
	  logf = - 0.5 * x
	         - 0.5 * (y[t].^2) * exp(-x)
			 - 0.5 * c2;
	}
  }
  if(isp == 0){
    adFunc[0] = logf;
  }
  else{
	adFunc[0] = logf - log(densn((x-mh)/vh));
  }
  return 1;
}

logph(const dht, const dh, const dphi, const dsig2, const dmu, const dt){
  decl f, c1, c2;
  if(dt == 0){
	c1 = dht|dh[1] - dmu;
	c2 = ((c1[1] - dphi*c1[0]).^2)/dsig2;
	f = - 0.5 * dht
	    - 0.5 * (y[dt].^2) * exp(-dht)
		- 0.5 * (1-dphi^2) / dsig2 * (c1[0].^2)
		- 0.5 * c2; 
  }
  else{
    if(dt != n-1){
	  c1 = dh[dt-1] | dht | dh[dt+1] - dmu;
	  c2 = ((c1[1:2] - dphi*c1[0:1]).^2)/dsig2;
	  f = - 0.5 * dht
	      - 0.5 * (y[dt].^2) * exp(-dht)
		  - 0.5 * sumc(c2);
	}
	else{
	  c1 = dh[dt-1] | dht - dmu;
	  c2 = ((c1[1] - dphi*c1[0]).^2)/dsig2;
	  f = - 0.5 * dht
	      - 0.5 * (y[dt].^2) * exp(-dht)
		  - 0.5 * c2;
	}
  }
  f = f - log(densn((dht-mh)/vh)) - logc;
  return(f);
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