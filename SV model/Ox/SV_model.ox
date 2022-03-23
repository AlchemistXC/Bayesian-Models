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
decl modephi, modeh, mlphi, mlh, Hphi, Hh, vphi2, vphi, vh2, vh, mphi, mh, sphi, sh;

progress_bar(const i, const iter, const burnin);
ranmu(const dphi, const dsig2, const dh);
ransigma2(const dphi, const dmu, const dh);
rantrancenorm(const r, const c, const mu, const sigma, const a, const b);
flphi(const dX, const adFunc, const avScore, const amHess);
flh(const dX, const adFunc, const avScore, const amHess);

main(){
  decl iter, burnin, irs, sample, sampleh;
  decl acphi, ach, car_count_h, car_count_phi;
  decl logrhophi, logrhoh, rhophi, rhoh, prop, logu, u, score, H, qphi1, qphi2, qh1, qh2, qphi, qh,
       pphi1, pphi2, ph1, ph2, pphi, ph;

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
  acphi = 0; // acceptance count of phi
  ach = 0;	// all acceptance count of h
  car_count_phi = 0;
  car_count_h = 0;
  for(decl i = -burnin; i < iter; i++){
	progress_bar(i, iter, burnin);
	// Simulate mu
	mu = ranmu(phi, sig2, h);
	// Simulate sigma2
	sig2 = ransigma2(phi, mu, h);
	// Simulate phi	(AR-MH)
	modephi = 0; // initial value
	MaxBFGS(flphi, &modephi, &mlphi, 0, 1);
	Num1Derivative(flphi, modephi, &sphi);
	Num2Derivative(flphi, modephi, &Hphi);
	vphi2 = -1/Hphi;
	vphi = sqrt(vphi2);
	mphi = modephi + vphi*vphi*sphi;
	// AR step
	do{
	  prop = rantrancenorm(1, 1, mphi, vphi, -1, 1);
	  flphi(prop, &qphi1, 0, 0);
	  qphi2 = mlphi + sphi * (prop - modephi) + 0.5 * Hphi * (prop - modephi)^2; 
	  qphi = qphi1 - qphi2; 
	  logu = log(ranu(1, 1));
	  car_count_phi++;	  
	}while(logu > qphi);
    //MH step
	flphi(phi, &pphi1, 0, 0);
	pphi2 = mlphi + sphi*(phi - modephi) + 0.5 * Hphi * (phi - modephi)^2;
	pphi = pphi1 - pphi2;
	if(pphi < 0){logrhophi = 0;}
	else{
	  if(qphi < 0){logrhophi = -pphi;}
	  else{logrhophi = qphi - pphi;}
	}
    logu = log(ranu(1, 1));
	if(logu <= logrhophi){
	  phi = prop;
	  acphi++;
	}
	//Simulate h (AR-MH)
    for(t = 0; t < n; t++){
	  MaxBFGS(flh, &modeh, &mlh, 0, 1);
	  Num1Derivative(flh, modeh, &sh);
	  Num2Derivative(flh, modeh, &Hh);
	  vh2 = -1/Hh;
	  vh = sqrt(vh2);
	  mh = modeh + vh*vh*sh;
	//A-R step
	  do{
		prop = mh + vh * rann(1, 1);
	    flh(prop, &qh1, 0, 0);
		qh2 = mlh + sh*(prop - modeh) + 0.5 * Hh * (prop - modeh)^2;
		qh = qh1 - qh2;
		logu = log(ranu(1, 1));
	  }while(logu > qh);
	//M-H step
	  flh(h[t], &ph1, 0, 0);
	  ph2 = mlh + sh*(h[t] - modeh) + 0.5 * Hh * (h[t] - modeh)^2;
	  ph = ph1 - ph2;
	  if(ph < 0){logrhoh = 0;}
	  else{
	    if(qh < 0){logrhoh = - ph;}
		else{logrhoh =qh - ph;}
	  }
	  logu = log(ranu(1, 1));
	  if(logu <= logrhoh){
	    h[t] = prop;
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

flphi(const dX, const adFunc, const avScore, const amHess){
  decl dphi, c1, c2, c3, c4, logf;
  dphi = dX;
  c1 = 1 - dphi^2;
  c2 = h - mu;
  c3 = sumc((c2[1:(n-1)] - dphi*c2[0:(n-2)]).^2); 
  logf =   log(densbeta((dphi+1)/2,a0,b0))
         + 0.5 * log(c1)
		 - 0.5 * c3 / sig2
		 - 0.5 * (c2[0].^2) * c1 / sig2;
  adFunc[0] = logf;
  return 1;			   
}

flh(const dX, const adFunc, const avScore, const amHess){
  decl dht, f, dmut, dsig2t, c;
  dht = dX;
  f = - y[t]^2 * exp(-dht) - 0.5 * dht;
  if(t == 0){
    dsig2t = sig2;
	dmut = (1 - phi) * mu + phi * h[1]; 
  }
  else{
    if(t != n-1){
	  c = 1 + phi^2;
	  dsig2t = sig2 / c;
	  dmut = phi * (h[t-1] + h[t+1])/ c;  
	}
	else{
	  dsig2t = sig2;
	  dmut = phi * h[n-2] + (1 - phi) * mu;
	}
  }
  adFunc[0] = f - 0.5 * (dht - dmut) * (dht - dmut) / dsig2t;
  return 1;
}