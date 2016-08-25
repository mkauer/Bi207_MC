/***********************************************************
A PRE-LOADED CONFIG FILE
------------------------------------------------------------

VERSION: 09.07.19

change log: (most recent at top)
-----------------------------------
19jul ~ 'ADC_offset' variable more accurately renamed
        to 'scint_quench'
29may + added xmove and ymove parameters
04nov + memory for some strings to use for runsDB
02nov + added some memory mappings to use later

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/

ofstream outfile;
Float_t calib;

const Int_t xsize=700,ysize=500;
const Int_t xmove=600,ymove=0;

const Int_t MAXPARS=25;

TF1 *func;
Double_t global_pars[MAXPARS];

Float_t N500=1,cedge500=1,csigma500=1,c2sigma500=1,
  Ngaus500=1,gausmean500=1,gaussig500=1,Nratio500=1,
  N1000=1,cedge1000=1,csigma1000=1,c2sigma1000=1,
  Ngaus1000=1,gausmean1000=1,gaussig1000=1,Nratio1000=1,
  Nbi500=1,cal482=1,Ek482=1,nEk482=1,s482=1,ns482=1,Ek482_err=1,
  s482_err=1,cal482_err=1,chi2_482=1,ndf_482=1,
  Nbi1000=1,cal976=1,Ek976=1,nEk976=1,s976=1,ns976=1,Ek976_err=1,
  s976_err=1,cal976_err=1,chi2_976=1,ndf_976=1,sys_s976=1,
  scint_quench=1,calib_err=1,tot_chi2=1,tot_ndf=1;

//TString rootfile;
string describe;
string rootfile;
Int_t daq,chan,maxbin,rebin;
Float_t gamma,beta,Eloss;
Bool_t lowres,thick;
char *pedrun,*betarun,*gammarun,*fitopt;
//char *fitopt;
//string pedrun,betarun,gammarun;

TFile *file;
TH1F *ptemp;
TH1F *gtemp;
TH1F *btemp;

TH1F *hped,*hgamma,*hbeta;

