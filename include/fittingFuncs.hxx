/***********************************************************
A BUNCH OF FUNCTIONS
------------------------------------------------------------

VERSION: 2012.10.25

change log: (most recent at top)
-------------------------------------------------
25oct12 + K,L,M shell independant fits
23sep09 ~ cleaned up a little
01sep09 + PMT gain power law
31aug09 + inverse expodential
04nov08 ~ this needs to get cleaned up soon

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/

// basic guassian
Double_t gg(Double_t x, Double_t mean, Double_t sigma)
{
  return TMath::Sqrt(1./(2.*TMath::Pi()))*(1./sigma)*TMath::Exp(-(x-mean)*(x-mean)/(2.*sigma*sigma));
}

// triple gaussian for 976kev
Double_t fitf(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return (par[0]*(7.03*gg(x,par[1],par[2])+1.84*gg(x,par[1]+72.3*par[3],par[2]*TMath::Sqrt(1+73*par[3]/par[1]))+0.545*gg(x,par[1]+84.3*par[3],par[2]*TMath::Sqrt(1+85*par[3]/par[1]))))+par[4]; 
}

//--- triple gaussian for 482kev
Double_t f_bi500(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return par[0]*(1.52*gg(x,par[1],par[2])+0.438*gg(x,par[1]+72.1*par[3],par[2]*TMath::Sqrt(1+73*par[3]/par[1]))+0.147*gg(x,par[1]+84.1*par[3],par[2]*TMath::Sqrt(1+85*par[3]/par[1])));
}

//--- triple gaussian for 976kev
Double_t f_bi1000(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return par[0]*(7.03*gg(x,par[1],par[2])+1.84*gg(x,par[1]+72.3*par[3],par[2]*TMath::Sqrt(1+73*par[3]/par[1]))+0.545*gg(x,par[1]+84.3*par[3],par[2]*TMath::Sqrt(1+85*par[3]/par[1])));
}

//--- full bi207 specturm
Double_t f_all(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return par[0]*TMath::Exp((-x)*par[1]) + f_compton(xval,&par[2]) + f_compton(xval,&par[6]) + f_bi500(xval,&par[10]) + f_bi1000(xval,&par[14]) + par[18];
}

//--- basic expodential
Double_t f_expo(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return (par[0]*TMath::Exp((-x)*par[1])) + par[2];
}

//--- basic inverse expodential
Double_t f_expo_inv(Double_t *xval, Double_t *par)
{
  return xval[0] + (par[0]*TMath::Exp((xval[0])*par[1])) - par[2];
}

//--- PMT gain power law
Double_t f_pmgain(Double_t *xval, Double_t *par)
{
  return (par[0]*TMath::Power(xval[0],par[1]));
}

//--- basic line
Double_t f_line(Double_t *xval, Double_t *par)
{
  return (par[1]*xval[0]) + par[0];
}

//--- basic guassian
Double_t f_gaus(Double_t *xval, Double_t *par) 
{
  Double_t x=xval[0];
  return par[0]*TMath::Sqrt(1/(2*TMath::Pi()))*(1/par[2])*TMath::Exp(-(x-par[1])*(x-par[1])/(2*par[2]*par[2]));
}

//--- basic compton
Double_t f_compton(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return par[0]*(1/(1+TMath::Exp((x-par[1])/par[2])))*(1+TMath::Exp(TMath::Power(x,1)*par[3]));
}

//--- double compton
Double_t f_gamma(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return (par[0]*TMath::Exp((-x*par[1]))) + f_compton(xval,&par[2]) + f_compton(xval,&par[6]) + par[10];
}

//--- MODIFIED double compton 2
Double_t f_gamma2(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return par[0]*TMath::Exp((-x*par[1])) + f_compton(xval,&par[2]) + f_compton(xval,&par[6]) + f_gaus(xval,&par[10]) + par[13];
}

//--- MODIFIED double compton + 1 gaussian
Double_t f_gamma3(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  Double_t Scale=par[6]*par[10];  //par[10] is now the ratio between compton and gaus normalization factors!
  return par[0]*TMath::Exp((-x*par[1])) + par[2]*(1/(1+TMath::Exp((x-par[3])/par[4])))*(1+TMath::Exp(TMath::Power(x,1)*par[5])) + par[6]*(1/(1+TMath::Exp((x-par[7])/par[8])))*(1+TMath::Exp(TMath::Power(x,1)*par[9])) + Scale*TMath::Sqrt(1/(2*TMath::Pi()))*(1/par[12])*TMath::Exp(-(x-par[11])*(x-par[11])/(2*par[12]*par[12])) + par[13];
}

//--- MODIFIED double compton + 2 gaussians
Double_t f_gamma4(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  Double_t Scale1=par[2]*par[6];
  Double_t Scale2=par[9]*par[13];  //par[10] is now the ratio between compton and gaus normalization factors!
  return par[0]*TMath::Exp((-x*par[1])) + par[2]*(1/(1+TMath::Exp((x-par[3])/par[4])))*(1+TMath::Exp(TMath::Power(x,1)*par[5])) + Scale1*gg(x,par[7],par[8]) + par[9]*(1/(1+TMath::Exp((x-par[10])/par[11])))*(1+TMath::Exp(TMath::Power(x,1)*par[12])) + Scale2*TMath::Sqrt(1/(2*TMath::Pi()))*(1/par[15])*TMath::Exp(-(x-par[14])*(x-par[14])/(2*par[15]*par[15])) + par[16];
}

//--- MODIFIED full bi207 specturm 2
Double_t f_all2(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  return par[0]*TMath::Exp((-x)*par[1]) + f_compton(xval,&par[2]) + f_compton(xval,&par[6]) + f_gaus(xval,&par[10]) + f_bi500(xval,&par[13]) + f_bi1000(xval,&par[17]) + par[21];
}

//--- MODIFIED full bi207 specturm + 1 compton gaussian
Double_t f_all3(Double_t *xval, Double_t *par)
{ 
  Double_t x=xval[0];
  Double_t Scale=par[6]*par[10];  //par[10] is now the ratio between compton and gaus normalization factors!
  return par[0]*TMath::Exp((-x*par[1])) + par[2]*(1/(1+TMath::Exp((x-par[3])/par[4])))*(1+TMath::Exp(TMath::Power(x,1)*par[5])) + par[6]*(1/(1+TMath::Exp((x-par[7])/par[8])))*(1+TMath::Exp(TMath::Power(x,1)*par[9])) + Scale*TMath::Sqrt(1/(2*TMath::Pi()))*(1/par[12])*TMath::Exp(-(x-par[11])*(x-par[11])/(2*par[12]*par[12])) + par[13]*(1.52*gg(x,par[14],par[15])+0.438*gg(x,par[14]+72.1*par[16],par[15]*TMath::Sqrt(1+73*par[16]/par[14]))+0.147*gg(x,par[14]+84.1*par[16],par[15]*TMath::Sqrt(1+85*par[16]/par[14]))) + par[17]*(7.03*gg(x,par[18],par[19])+1.84*gg(x,par[18]+72.3*par[20],par[19]*TMath::Sqrt(1+73*par[20]/par[18]))+0.545*gg(x,par[18]+84.3*par[20],par[19]*TMath::Sqrt(1+85*par[20]/par[18]))) + par[21];
}

//--- MODIFIED full bi207 specturm + 2 compton gaussians
Double_t f_all4(Double_t *xval, Double_t *par)
{ 
  Double_t x=xval[0];
  
  Double_t exp0=par[0];
  Double_t exp1=par[1];
  
  Double_t Cn500=par[2];
  Double_t Ce500=par[3];
  Double_t Cs500=par[4];
  Double_t Css500=par[5];
  Double_t CGn500=par[6];
  Double_t CGe500=par[7];
  Double_t CGs500=par[8];
  
  Double_t Cn1000=par[9];
  Double_t Ce1000=par[10];
  Double_t Cs1000=par[11];
  Double_t Css1000=par[12];
  Double_t CGn1000=par[13];
  Double_t CGe1000=par[14];
  Double_t CGs1000=par[15];
  
  Double_t Gn482=par[16];
  Double_t Ge482=par[17];
  Double_t Gs482=par[18];
  Double_t Gc482=par[19];
  
  Double_t Gn976=par[20];
  Double_t Ge976=par[21];
  Double_t Gs976=par[22];
  Double_t Gc976=par[23];
  
  Double_t offset=par[24];
  Double_t Scale1=Cn500*CGn500;
  Double_t Scale2=Cn1000*CGn1000;
  
  Double_t Exp=0,Comp500=0,CompG500=0,Comp1000=0,CompG1000=0;
  Double_t Bi482=0,Bi482k=0,Bi482l=0,Bi482m=0;
  Double_t Bi976=0,Bi976k=0,Bi976l=0,Bi976m=0;
  
  Exp = exp0*TMath::Exp((-x*exp1));
  
  Comp500   = Cn500*(1/(1+TMath::Exp((x-Ce500)/Cs500)))*(1+TMath::Exp(TMath::Power(x,1)*Css500));
  CompG500  = Scale1*TMath::Gaus(x,CGe500,CGs500);
  Comp1000  = Cn1000*(1/(1+TMath::Exp((x-Ce1000)/Cs1000)))*(1+TMath::Exp(TMath::Power(x,1)*Css1000));
  CompG1000 = Scale2*TMath::Gaus(x,CGe1000,CGs1000);
  
  Bi482k = 1.520*TMath::Gaus(x,Ge482,Gs482);
  Bi482l = 0.438*TMath::Gaus(x,Ge482+72.1*Gc482,Gs482*TMath::Sqrt(1+73*Gc482/Ge482));
  Bi482m = 0.147*TMath::Gaus(x,Ge482+84.1*Gc482,Gs482*TMath::Sqrt(1+85*Gc482/Ge482));
  Bi482  = Gn482*(Bi482k + Bi482l + Bi482m);
  
  Bi976k = 7.030*TMath::Gaus(x,Ge976,Gs976);
  Bi976l = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  Bi976m = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  Bi976  = Gn976*(Bi976k + Bi976l + Bi976m);

  return Exp + Comp500 + CompG500 + Bi482 + Comp1000 + CompG1000 + Bi976 + offset;
}

//--- Delta-E Fit to the Bi207
Double_t deltaE(Double_t *xval, Double_t *par)
{ 
  Double_t x=xval[0];
  return par[0]*(1.52*gg(x,par[1],par[2])+0.438*gg(x,par[1]+72.1*par[3],par[2]*TMath::Sqrt(1+73*par[3]/par[1]))+0.147*gg(x,par[1]+84.1*par[3],par[2]*TMath::Sqrt(1+85*par[3]/par[1]))) + par[4]*(7.03*gg(x,par[5],par[6])+1.84*gg(x,par[5]+72.3*par[7],par[6]*TMath::Sqrt(1+73*par[7]/par[5]))+0.545*gg(x,par[5]+84.3*par[7],par[6]*TMath::Sqrt(1+85*par[7]/par[5]))) + par[8];
}

//--- Fit to the Bi207 976keV with Landau tail
Double_t f_landau_976(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  Double_t Ln976=par[4];
  Double_t El976=par[5];
  Double_t Ls976=par[6];
  
  Double_t offset=par[7];
  
  Double_t mpv976=0,land976=0;
  Double_t Bi976=0,Bi976k=0,Bi976l=0,Bi976m=0;
  
  mpv976  = Ge976-(El976*Gc976);
  land976 = Ln976*(TMath::Landau((2*mpv976)-x,mpv976,Ls976));
  
  Bi976k  = 7.030*TMath::Gaus(x,Ge976,Gs976);
  Bi976l  = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  // testing
  //Bi976l  = par[8]*TMath::Gaus(x,Ge976+par[9]*Gc976,Gs976*TMath::Sqrt(1+par[9]*Gc976/Ge976));
  Bi976m  = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  Bi976   = Gn976*(Bi976k + Bi976l + Bi976m);
    
  return land976 + Bi976 + offset;
}

//--- Fit to the Bi207 976keV and 482keV with Landau tails
Double_t f_landau_full(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  Double_t Ln976=par[4];
  Double_t El976=par[5];
  Double_t Ls976=par[6];
  
  Double_t offset=par[7];
  
  Double_t Gn482=par[8];
  Double_t Ge482=par[9];
  Double_t Gs482=par[10];
  Double_t Gc482=par[11];
  Double_t Ln482=par[12];
  Double_t El482=par[13];
  Double_t Ls482=par[14];
  
  Double_t mpv976=0,land976=0,gaus976=0;
  Double_t mpv482=0,land482=0,gaus482=0;
  Double_t Bi482=0,Bi482k=0,Bi482l=0,Bi482m=0;
  Double_t Bi976=0,Bi976k=0,Bi976l=0,Bi976m=0;
  
  mpv482  = Ge482-(El482*Gc482);
  land482 = Ln482*(TMath::Landau((2*mpv482)-x,mpv482,Ls482));
  Bi482k  = 1.520*TMath::Gaus(x,Ge482,Gs482);
  Bi482l  = 0.438*TMath::Gaus(x,Ge482+72.1*Gc482,Gs482*TMath::Sqrt(1+73*Gc482/Ge482));
  Bi482m  = 0.147*TMath::Gaus(x,Ge482+84.1*Gc482,Gs482*TMath::Sqrt(1+85*Gc482/Ge482));
  Bi482   = Gn482*(Bi482k + Bi482l + Bi482m);
    
  mpv976  = Ge976-(El976*Gc976);
  land976 = Ln976*(TMath::Landau((2*mpv976)-x,mpv976,Ls976));
  Bi976k  = 7.030*TMath::Gaus(x,Ge976,Gs976);
  Bi976l  = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  Bi976m  = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  Bi976   = Gn976*(Bi976k + Bi976l + Bi976m);
    
  return Bi976+land976+offset+Bi482+land482;
}

// POISSON
Double_t f_poisson(Double_t *xval, Double_t *par)
{
  return (par[0]*TMath::Poisson(xval[0]/par[2],par[1]))+(par[3]*TMath::Gaus(xval[0],0.0,par[4]));
}

// POISSON 2
Double_t f_poisson2(Double_t *xval, Double_t *par)
{
  return (par[0]*TMath::Poisson(xval[0]/par[2],par[1]))+(par[3]*TMath::Gaus(xval[0]/par[2],0.0,par[4]));
}

// 976 keV fit 
Double_t f_976kev(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  
  Double_t Bi976=0,Bi976k=0,Bi976l=0,Bi976m=0;
  
  Bi976k  = 7.030*TMath::Gaus(x,Ge976,Gs976);
  Bi976l  = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  Bi976m  = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  Bi976   = Gn976*(Bi976k + Bi976l + Bi976m);
    
  return Bi976;
}
// 976 K-shell 
Double_t f_976K(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  
  Double_t Bi976k=0;
  
  Bi976k  = 7.030*TMath::Gaus(x,Ge976,Gs976);
  //Bi976l  = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  //Bi976m  = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  //Bi976   = Gn976*(Bi976k + Bi976l + Bi976m);
    
  return Gn976*Bi976k;
}
// 976 L-shell 
Double_t f_976L(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  
  Double_t Bi976l=0;
  
  //Bi976k  = 7.030*TMath::Gaus(x,Ge976,Gs976);
  Bi976l  = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  //Bi976m  = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  //Bi976   = Gn976*(Bi976k + Bi976l + Bi976m);
    
  return Gn976*Bi976l;
}
// 976 M-shell 
Double_t f_976M(Double_t *xval, Double_t *par)
{
  Double_t x=xval[0];
  
  Double_t Gn976=par[0];
  Double_t Ge976=par[1];
  Double_t Gs976=par[2];
  Double_t Gc976=par[3];
  
  Double_t Bi976m=0;
  
  //Bi976k  = 7.030*TMath::Gaus(x,Ge976,Gs976);
  //Bi976l  = 1.840*TMath::Gaus(x,Ge976+72.3*Gc976,Gs976*TMath::Sqrt(1+73*Gc976/Ge976));
  Bi976m  = 0.545*TMath::Gaus(x,Ge976+84.3*Gc976,Gs976*TMath::Sqrt(1+85*Gc976/Ge976));
  //Bi976   = Gn976*(Bi976k + Bi976l + Bi976m);
    
  return Gn976*Bi976m;
}


//--- my sweet ass SPE fit
Double_t f_spe(Double_t *xval, Double_t *par)
{
  const Double_t sr2pi=TMath::Sqrt(2*TMath::Pi());
  
  Double_t x=xval[0];
  Double_t N=par[0];
  Double_t ped_mean=par[1];
  Double_t ped_sig=par[2];
  Double_t Npe=par[3];
  Double_t pe_mean=par[4];
  Double_t pe_sig=par[5];
  Double_t Df=par[6];
  Double_t Ds=par[7];
    
  Double_t ped_term,pe_term,dyn_term,sum;
  
  ped_term = N*TMath::Exp(-Npe)
    *(1./(sr2pi*ped_sig))
    *TMath::Exp(-(TMath::Power(x-ped_mean,2))/(2.*ped_sig*ped_sig));
  
  sum=0;
  for(Int_t n=1; n<=12; n++){
    sum += TMath::Exp(-Npe)
      *TMath::Power(Npe,n)
      *TMath::Exp(-TMath::Power(x-ped_mean-Double_t(n)*pe_mean,2)/(2.*(ped_sig*ped_sig+Double_t(n)*pe_sig*pe_sig)))
      *(1./TMath::Factorial(n))
      *(1./(sr2pi*TMath::Sqrt(2.*(ped_sig*ped_sig+Double_t(n)*pe_sig*pe_sig))));
  }
  
  pe_term = N*(1-Df)*sum;
  
  dyn_term = N*Df*(1-TMath::Exp(-Npe))
    *TMath::Exp(-(TMath::Power((x-ped_mean-pe_mean)/Ds,2))/(2.*((ped_sig*ped_sig+pe_sig*pe_sig)/(Ds*Ds))))
    *(1./(sr2pi*TMath::Sqrt((ped_sig*ped_sig+pe_sig*pe_sig)/(Ds*Ds))));
  
  return ped_term + pe_term + dyn_term;

}

