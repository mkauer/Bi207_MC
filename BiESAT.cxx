/***********************************************************
BiESAT = Bi207 ENERGY SPECTRA ANALYSIS TOOLKIT
------------------------------------------------------------
This generally works for all bi207 spectra with exceptions:
- This cannot handle changing the x-axis into units of MeV
  but units of keV are fine (the applies more to fitting of
  simulated data).

VERSION: 2012.10.25

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/
gROOT->Reset();

#include "./include/fittingFuncs.hxx"
#include "./include/init_memory.cxx"
#include "./include/manip_files.hxx"
#include "./include/manip_params.hxx"
#include "./include/manip_histos.hxx"

#include "./quickfit_976kev.cxx"
#include "./quickfit_landau_976.cxx"

Int_t verbose=1;
Double_t eg_diff=576; // diff between 500keV compton edge and the 976keV electron.

Int_t BiESAT(string runnumber,Int_t DB=0,string breakpoint=""){
  string workdir="";
  if(DB==0) workdir="./DBruns/";
  if(DB==1) workdir="./DBsimu/";
  if(DB==2) workdir="./DBtemp/";
  string configfile=runnumber+".conf";
  rootfile=runnumber+".root";
  
  cout<<"Looking in ==> "<<workdir<<" for config file ==> "<<configfile<<endl;
  gROOT->LoadMacro((workdir+configfile).c_str());
  config();
  
  file=new TFile((workdir+rootfile).c_str(),"READ");
  ptemp=(TH1F*)file->Get(pedrun);
  gtemp=(TH1F*)file->Get(gammarun);
  btemp=(TH1F*)file->Get(betarun);
  
  outfile.open("AA_RESULTS.txt");
  
  beautify();
  
  cout<<"\n\n ===>  "<<describe<<endl;
  cout<<"\n\t+++++++++++++++++++++++++++++++++++++++++++ \n";
  if(lowres)  cout<<"\t  USING FIT MODE  ==>  low-resolution \n";
  if(!lowres) cout<<"\t  USING FIT MODE  ==>  hi-resolution \n";
  if(thick)   cout<<"\t  USING FIT MODE  ==>  thick scint \n";
  if(!thick)  cout<<"\t  USING FIT MODE  ==>  thin scint \n";
  cout<<"\t+++++++++++++++++++++++++++++++++++++++++++ \n";
  
  calib=(beta-gamma)/(eg_diff-Eloss);
  if(calib<=0){
    cout<<"\n\n\t CALIB IS NEGATIVE! ASSUMING ABSOLUTE VALUE! \n\n\n";
    calib=TMath::Abs(calib);
  }
  cout<<"\t  Preliminary Calib    ==>  "<<calib<<"\n";
  cout<<"\t  Assumed Energy Loss  ==>  "<<Eloss<<"\n";
  cout<<"\t+++++++++++++++++++++++++++++++++++++++++++ \n\n\n";
  
  Int_t xmin=0,xmax=maxbin,pmin=0,pmax=maxbin;
  TCanvas *ctmp = new TCanvas("ctmp","ctmp",1,1,10,10);
  ctmp->cd();
  
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //   FIT THE PEDISTAL WITH BASIC GAUSSIAN
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  
  Float_t ped=0;
  Float_t sped=0;
  Float_t ped_err=0;
  Float_t sped_err=0;
  
  if(ptemp){
    hped=(TH1F*)resizeHist(ptemp,0,maxbin,"hped",pedrun);
    
    TF1 *pedfit = new TF1("pedfit","gaus",xmin,xmax);
    pedfit->SetLineColor(2);
    pedfit->SetLineWidth(2);
    pedfit->SetParLimits(0,0,1e7);
    pedfit->SetParLimits(1,0,1e3);
    pedfit->SetParLimits(2,0,1e2);
        
    TCanvas *c70 = new TCanvas("c70","- Pedistal Fit -",0,0,xsize,ysize);
    c70->cd();
    hped->Draw();
    hped->SetAxisRange(2,1000,"X");
    ped=hped->GetMaximumBin();
    sped=hped->GetRMS();
    xmin=0;
    xmax=ped+(1*sped);
    hped->Fit("pedfit","Q","",xmin,xmax);
    
    ped = pedfit->GetParameter(1);
    sped = pedfit->GetParameter(2);
    xmin=ped-(3*sped);
    xmax=ped+(2*sped);
    hped->Fit("pedfit","QWW","",xmin,xmax);
    
    ped = pedfit->GetParameter(1);
    sped = pedfit->GetParameter(2);
    ped_err = pedfit->GetParError(1);
    sped_err = pedfit->GetParError(2);
    
    cout<<"\n\n\t Pedestal Mean = "<<ped<<"  /  Sigma = "<<sped<<"\n";
    
    pmin=ped-8*sped;
    pmax=ped+10*sped;
    if(pmin<0) pmin=0;
    hped->SetAxisRange(pmin,pmax,"X");
    c70->Update();
    c70->Print("AA_pedistal.png");
    
    ctmp->cd();
    delete pedfit;
  }
  
  const Int_t hmin = 1;
  const Int_t hmax = maxbin-ped-1;
  
  if(!ptemp){ 
    cout<<"\n\t NO PEDESTAL FILE FOUND ==> using a vaule = "<<ped<<"\n\n";
  }
  
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //   CONVERTING, SHIFTING, REBINNING AND NORMALIZING
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  
  // CONVERTING
  shiftHist(btemp,ped);
  TH1F *hbeta=(TH1F*)resizeHist(btemp,0,maxbin,"hbeta",betarun);
  hbeta->SetXTitle("QDC bins");
  hbeta->SetYTitle("counts/bin");
  hbeta->SetTitle(describe.c_str());
  
  shiftHist(gtemp,ped);
  TH1F *hgamma=(TH1F*)resizeHist(gtemp,0,maxbin,"hgamma",gammarun);
  hgamma->SetXTitle("QDC bins");
  hgamma->SetYTitle("counts/bin");
  hgamma->SetTitle(describe.c_str());
  
  TF1 *cfit = new TF1("cfit","gaus",0,maxbin);
  cfit->SetLineColor(2);
  cfit->SetLineWidth(2);
  
  // SHIFTING
  Int_t newrebin=rebin+4;
  if(verbose) cout<<"\n\n\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  btemp->Rebin(newrebin);
  xmin=340*calib;
  xmax=440*calib;
  btemp->SetAxisRange(xmin,xmax,"X");
  if(verbose) cout<<"\t looking for BETA max bin between "<<xmin<<" - "<<xmax<<" bins \n";
  Double_t comp_beta=btemp->GetMaximumBin()*(newrebin);
  btemp->SetAxisRange(0,maxbin,"X");
  xmin=comp_beta-30*calib;
  xmax=comp_beta+30*calib;
  if(verbose) cout<<"\t looking for BETA 500keV compton between "<<xmin<<" - "<<xmax<<" bins \n";
  btemp->Fit("cfit","Q","",xmin,xmax);
  comp_beta=cfit->GetParameter(1);
  
  gtemp->Rebin(newrebin);
  xmin=340*calib;
  xmax=440*calib;
  gtemp->SetAxisRange(xmin,xmax,"X");
  if(verbose) cout<<"\t looking for GAMMA max bin between "<<xmin<<" - "<<xmax<<" bins \n";
  Double_t comp_gamma=gtemp->GetMaximumBin()*(newrebin);
  gtemp->SetAxisRange(0,maxbin,"X");
  xmin=comp_gamma-30*calib;
  xmax=comp_gamma+30*calib;
  if(verbose) cout<<"\t looking for GAMMA 500keV compton between "<<xmin<<" - "<<xmax<<" bins \n";
  gtemp->Fit("cfit","Q","",xmin,xmax);
  comp_gamma=cfit->GetParameter(1);
  if(verbose) cout<<"\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
  
  cout<<"\n\t SHIFTING --> "<<gammarun<<" --> "<<comp_gamma-comp_beta<<" bins\n\n";
  shiftHist(hgamma,comp_gamma-comp_beta);
  
  if(verbose){
    TCanvas *look = new TCanvas("look","look",0,0,xsize,ysize);
    look->Divide(1,2,0.02,0.02,0);
    look->cd(1);
    btemp->Draw();
    look->Update();
    look->cd(2);
    gtemp->Draw();
    look->Update();
    ctmp->cd();
  }
  delete cfit;
  //break;
  
  // REBINNING
  hgamma->Rebin(rebin);
  hbeta->Rebin(rebin);
  
  // NORMALIZING
  Int_t normv=0;   /// if normalization fails, try setting this to 1
  normalizer(hgamma,hbeta,normv);
  
  drawMix(hgamma,hbeta);
  ctmp->cd();
  drawSub(hgamma,hbeta,2);
  ctmp->cd();
  TH1F *sub=(TH1F*)subtract(hgamma,hbeta);
  
  if(breakpoint=="sub"){
    cout<<"\n\nBiESAT OPTION: stopping after histogram subtraction \n\n\n";
    delete ctmp;
    break;
  }
  
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //     FIT THE COMPTON PEAKS 
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  
  delete func;
  func = new TF1("func",f_all4,hmin,hmax,MAXPARS);
  setParNames(func,"bi207");
  func->SetLineWidth(2);
  func->SetLineColor(2);
  
  if(verbose) cout<<"\n\n\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
  
  //---- expo
  func->SetParameter(0,1.0e2);
  func->SetParLimits(0,1.0,1.0e5);
  func->SetParameter(1,1.0e-5);
  func->SetParLimits(1,0.0,1.0);
    
  //---- 500 compton
  func->SetParameter(2,200.1);
  func->SetParLimits(2,1.0,1.0e3);
  func->SetParameter(3,400.*calib);
  func->SetParLimits(3,350.*calib,420.*calib);
  if(verbose) cout<<"\t looking for 500keV compton between "<<360.*calib<<" - "<<420.*calib<<" bins \n";
  func->SetParameter(4,1.1);
  func->SetParLimits(4,1.0,100);
  func->SetParameter(5,0.0001);
  func->SetParLimits(5,0.0,0.02);
    
  //---- 500 multi-compton
  func->SetParameter(6,10.1);
  func->SetParLimits(6,1.0,1.0e4);
  func->SetParameter(7,440.*calib);
  func->SetParLimits(7,400.*calib,540.*calib); //--- SENSITIVE PARAMETER!!
  if(verbose && !lowres && thick) cout<<"\t looking for 500keV multi-compton between "<<400.*calib<<" - "<<540.*calib<<" bins \n";
  func->SetParameter(8,10.1);
  func->SetParLimits(8,5.0,100.0);
      
  //---- 1000 compton
  func->SetParameter(9,20.1);
  func->SetParLimits(9,1.0,1.0e3);
  func->SetParameter(10,900.*calib);
  func->SetParLimits(10,700.*calib,1100.*calib);
  if(verbose) cout<<"\t looking for 1000keV compton between "<<700.*calib<<" - "<<1100.*calib<<" bins \n";
  func->SetParameter(11,10.1);
  func->SetParLimits(11,1.0,100.0);
  func->SetParameter(12,0.00001);
  func->SetParLimits(12,0.0,0.01);
    
  //---- 1000 multi-compton
  func->SetParameter(13,10.1);
  func->SetParLimits(13,1.0,1.0e4);
  func->SetParameter(14,820.*calib);
  func->SetParLimits(14,780.*calib,850.*calib); //--- SENSITIVE PARAMETER!!
  if(verbose && !lowres) cout<<"\t looking for 1000keV multi-compton between "<<780.*calib<<" - "<<850.*calib<<" bins \n";
  func->SetParameter(15,20.1);
  func->SetParLimits(15,10.0,500.0);
  
  //---- 482 bi207
  func->FixParameter(16,0.0);
  func->FixParameter(17,1.0e-2);
  func->FixParameter(18,1.0e-2);
  func->FixParameter(19,1.0e-2);
  
  //---- 976 bi207
  func->FixParameter(20,0.0);
  func->FixParameter(21,1.0e-2);
  func->FixParameter(22,1.0e-2);
  func->FixParameter(23,1.0e-2);
  
  //---- offset
  func->SetParameter(24,1.1);
  func->SetParLimits(24,0.0,500.0);
  
  if(lowres || !thick){
    //---- 500 multi-compton
    func->FixParameter(6,0.0);
    func->FixParameter(7,1.0e-2);
    func->FixParameter(8,1.0e-2);
  }
    
  if(lowres){
    //---- 1000 multi-compton
    func->FixParameter(13,0.0);
    func->FixParameter(14,1.0e-2);
    func->FixParameter(15,1.0e-2);
  }
  
  if(verbose) cout<<"\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n";
  
  ///////////////////////////////////////////////////////////////////
  TCanvas *c71 = new TCanvas("c71","- Gamma Fit -",0,0,xsize,ysize);
  c71->cd();
  xmin=270.*calib;
  xmax=1200.*calib;
  hgamma->Draw();
  if(lowres || !thick) hgamma->Fit("func",fitopt,"",xmin,xmax);
  ///////////////////////////////////////////////////////////////////
  Bool_t failed=false;
  if(!lowres && thick){
    Bool_t goodfit=false;
    Int_t inc=0;
    while(!goodfit){
      cout<<"\n\n\t==============================================================\n";
      cout<<"\t Fitting the Gamma spectra for the "<<inc+1<<"-th time"<<endl;
      func->SetParameter(10,(900.-(inc*10.))*calib);
      func->SetParLimits(10,(780.-(inc*10.))*calib,(920.-(inc*10.))*calib);
      cout<<"\t Looking for the gamma ==> "<<(780.+(inc*10.))*calib<<" - "
	  <<(920.+(inc*10.))*calib<<endl;
      func->SetParameter(14,(900.+(inc*10.))*calib);
      func->SetParLimits(14,(860.+(inc*10.))*calib,(1020.+(inc*10.))*calib);
      cout<<"\t Looking for the multi-comp bump ==> "<<(860.+(inc*10.))*calib
	  <<" - "<<(1020.+(inc*10.))*calib<<endl;
      cout<<"\t==============================================================\n\n\n"; 
      hgamma->Fit("func",fitopt,"",xmin,xmax);
      if(func->GetParameter(10)>func->GetParameter(14)){
	goodfit=false;
	inc++;
      } 
      if(func->GetParameter(10)<func->GetParameter(14)){
	goodfit=true;
      } 
      if(inc>3){
	failed=true;
	break;
      }
    }
  }
  if(failed){
    cout<<"\n\n\t COULD NOT FIT GAMMA SPECTRA CORRECTLY!!! \n\n"<<endl;
    func->SetParLimits(10,700.*calib,1100.*calib);
    func->FixParameter(13,0.0);
    func->FixParameter(14,1.0e-2);
    func->FixParameter(15,1.0e-2);
    hgamma->Fit("func",fitopt,"",xmin,xmax);
  }
  ///////////////////////////////////////////////////////////////////
  
  N500         = func->GetParameter(2);
  cedge500     = func->GetParameter(3);
  csigma500    = func->GetParameter(4);
  c2sigma500   = func->GetParameter(5);
  Nratio500    = func->GetParameter(6);
  gausmean500  = func->GetParameter(7);
  gaussig500   = func->GetParameter(8);
    
  N1000        = func->GetParameter(9);
  cedge1000    = func->GetParameter(10);
  csigma1000   = func->GetParameter(11);
  c2sigma1000  = func->GetParameter(12);
  Nratio1000   = func->GetParameter(13);
  gausmean1000 = func->GetParameter(14);
  gaussig1000  = func->GetParameter(15);
  
  calib = (((cedge1000+gausmean1000)/2)-cedge500)/470.;
  if(thick) calib = (gausmean1000-cedge500)/510.;
  if(lowres) calib = (cedge1000-cedge500)/500.;
  
  Ek482        = (482.-Eloss)*calib;
  Ek976        = (976.-Eloss)*calib;
  s482         = (0.15*Ek482)/2.35;
  s976         = (0.15*Ek976)/2.35;
  
  tot_chi2=func->GetChisquare();
  tot_ndf=func->GetNDF();
  
  cout<<"\n\n++++++++++ COMPTON FIT RESULTS +++++++++++++\n";
  cout<<"\t N500         = "<<N500<<endl;
  cout<<"\t N1000        = "<<N1000<<endl;
  cout<<"\t cedge500     = "<<cedge500<<endl;
  cout<<"\t cedge1000    = "<<cedge1000<<endl;
  cout<<"\t csigma500    = "<<csigma500<<endl;
  cout<<"\t csigma1000   = "<<csigma1000<<endl;
  cout<<"\t c2sigma500   = "<<c2sigma500<<endl;
  cout<<"\t c2sigma1000  = "<<c2sigma1000<<endl;
  cout<<"\t Nratio500    = "<<Nratio500<<endl;
  cout<<"\t gausmean500  = "<<gausmean500<<endl;
  cout<<"\t gaussig500   = "<<gaussig500<<endl;
  cout<<"\t Nratio1000   = "<<Nratio1000<<endl;
  cout<<"\t gausmean1000 = "<<gausmean1000<<endl;
  cout<<"\t gaussig1000  = "<<gaussig1000<<endl;
  cout<<"\t Calib        = "<<calib;      
    
  cout<<"\n\n++++++++++ ESTIMATED PARAMETERS ++++++++++++\n";
  cout<<"\t Ek976        = "<<Ek976<<endl;
  cout<<"\t s976         = "<<s976<<endl;
  cout<<"\t Ek482        = "<<Ek482<<endl;
  cout<<"\t s482         = "<<s482<<"\n\n\n";
    
  pmin=xmin-3*csigma500;
  pmax=xmax+3*csigma1000;
  if(pmin<hmin) pmin=hmin;
  if(pmax>hmax) pmax=hmax;
  hgamma->SetAxisRange(pmin,pmax,"X");
  func->Draw("same");
  c71->Update();
  c71->Print("AA_gammas.png");
  c71->Print("AA_gammas.pdf");
  ctmp->cd();
  
  if(breakpoint=="gamma"){
    cout<<"\n\nBiESAT OPTION: stopping after fit to gamma spectrum \n\n\n";
    delete func;
    delete ctmp;
    break;
  }
  
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //      1st ROUGH FIT TO BI207 SPECTRUM 
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
    
  TH1F *hbeta1 = (TH1F*)cloneHist(hbeta,"Rough Fit 2",betarun);
  TH1F *hbeta2 = (TH1F*)cloneHist(hbeta,"Final Fit",betarun);
    
  //---- release parameters
  func->ReleaseParameter(16);
  func->ReleaseParameter(17);
  func->ReleaseParameter(18);
  func->ReleaseParameter(19);
  func->ReleaseParameter(20);
  func->ReleaseParameter(21);
  func->ReleaseParameter(22);
  func->ReleaseParameter(23);
    
  //---- expo
  //func->FixParameter(0,0.0);
  func->SetParameter(0,func->GetParameter(0));
  func->SetParLimits(0,0.0,1.0e5);
  func->FixParameter(1,func->GetParameter(1));
  
  //---- 500 compton
  func->SetParameter(2,N500);
  func->SetParLimits(2,0.1,1.0e4);
  func->FixParameter(3,cedge500);
  func->FixParameter(4,csigma500);
  func->FixParameter(5,c2sigma500);
  
  //---- 500 multi-compton
  func->SetParameter(6,Nratio500);
  func->SetParLimits(6,0.8*Nratio500,1.2*Nratio500);
  func->FixParameter(7,gausmean500);
  func->FixParameter(8,gaussig500);
  
  //---- 1000 compton
  func->FixParameter(9,N1000);
  func->FixParameter(10,cedge1000);
  func->FixParameter(11,csigma1000);
  func->FixParameter(12,c2sigma1000);
  
  //---- 1000 multi-compton
  func->FixParameter(13,Nratio1000);
  func->FixParameter(14,gausmean1000);
  func->FixParameter(15,gaussig1000);
  
  //---- 482 bi207
  func->SetParameter(16,1.0);
  func->SetParLimits(16,1.0,1.0e6);
  func->SetParameter(17,(482.-Eloss)*calib);
  func->SetParLimits(17,(440.-Eloss)*calib,(500.-Eloss)*calib);
  func->SetParameter(18,s482);
  func->SetParLimits(18,s482*0.20,s482*3.00);
  func->SetParameter(19,calib);
  func->SetParLimits(19,calib*0.2,calib*1.8);
    
  //---- 976 bi207
  func->SetParameter(20,1.0);
  func->SetParLimits(20,1.0,1.0e7);
  func->SetParameter(21,(976.-Eloss)*calib);
  func->SetParLimits(21,(900.-Eloss)*calib,(1000.-Eloss)*calib);  // bounds give troubles (Eloss problem???)
  func->SetParameter(22,s976);
  func->SetParLimits(22,s976*0.20,s976*3.0);
  func->SetParameter(23,calib);
  func->SetParLimits(23,calib*0.2,calib*1.8);
  
  //---- offset
  func->SetParameter(24,func->GetParameter(24));
  func->SetParLimits(24,0.0,500.0);
  
  if(lowres || !thick){
    //---- 500 multi-compton
    func->FixParameter(6,0.0);
    func->FixParameter(7,1.0e-2);
    func->FixParameter(8,1.0e-2);
  }
  
  if(lowres){
    //---- 1000 multi-compton
    func->FixParameter(13,0.0);
    func->FixParameter(14,1.0e-2);
    func->FixParameter(15,1.0e-2);
  }
  
  ///////////////////////////////////////////////////////////////////
  TCanvas *c72 = new TCanvas("c72","- Rough Fit 1 -",xmove,0,xsize,ysize);
  c72->cd();
  xmin=280*calib;
  xmax=1200*calib;
  if(xmax>hmax) xmax=hmax;
  hbeta->Draw();
  hbeta->Fit("func",fitopt,"",xmin,xmax);
  ///////////////////////////////////////////////////////////////////
  
  N500         = func->GetParameter(2);
  cedge500     = func->GetParameter(3);
  csigma500    = func->GetParameter(4);
  c2sigma500   = func->GetParameter(5);
  Nratio500    = func->GetParameter(6);
  gausmean500  = func->GetParameter(7);
  gaussig500   = func->GetParameter(8);
    
  N1000        = func->GetParameter(9);
  cedge1000    = func->GetParameter(10);
  csigma1000   = func->GetParameter(11);
  c2sigma1000  = func->GetParameter(12);
  Nratio1000   = func->GetParameter(13);
  gausmean1000 = func->GetParameter(14);
  gaussig1000  = func->GetParameter(15);
    
  Ek482        = func->GetParameter(17);
  Ek482_err    = func->GetParError(17);
  s482         = func->GetParameter(18);
  s482_err     = func->GetParError(18);
  cal482       = func->GetParameter(19);
  cal482_err   = func->GetParError(19);
  
  Ek976        = func->GetParameter(21);
  Ek976_err    = func->GetParError(21);
  s976         = func->GetParameter(22);
  s976_err     = func->GetParError(22);
  cal976       = func->GetParameter(23);
  cal976_err   = func->GetParError(23);
    
  calib        = (Ek976-Ek482)/494.;
  if(lowres || thick){
    if(lowres)calib = cal976;
    if(thick) calib = (Ek976-cedge500)/(eg_diff-Eloss);
    cal482     = calib;
    cal976     = calib;
  }
  
  scint_quench = 0;
  nEk976       = Ek976;
  nEk482       = Ek482;
  if(scint_quench>0){
    nEk976     = Ek976+scint_quench;
    nEk482     = Ek482+scint_quench;
  }
  ns976        = s976;
  if(sped<s976) ns976 = TMath::Sqrt((s976*s976)-(sped*sped));
  ns482        = s482;
  if(sped<s482) ns482 = TMath::Sqrt((s482*s482)-(sped*sped));
  
  tot_chi2=func->GetChisquare();
  tot_ndf=func->GetNDF();
    
  pmin=xmin-(3*csigma500);
  pmax=xmax+(4*s976);
  if(pmin<hmin) pmin=hmin;
  if(pmax>hmax) pmax=hmax;
  hbeta->SetAxisRange(pmin,pmax,"X");
  c72->Update();
  c72->Print("AA_bi207_rough1.png");
  ctmp->cd();
  
  cout<<"\n\n ROUGH FIT 1"<<endl;
  megaPrint();
  
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //      2nd ROUGH FIT TO BI207 SPECTRUM 
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
    
  //---- expo
  func->FixParameter(0,func->GetParameter(0));
  func->FixParameter(1,func->GetParameter(1));
  
  //---- 500 compton
  func->FixParameter(2,N500);
  func->FixParameter(3,cedge500);
  func->FixParameter(4,csigma500);
  func->FixParameter(5,c2sigma500);
    
  //---- 500 multi-compton
  func->FixParameter(6,Nratio500);
  func->FixParameter(7,gausmean500);
  func->FixParameter(8,gaussig500);
  
  //---- 1000 compton
  func->FixParameter(9,N1000);
  func->FixParameter(10,cedge1000);
  func->FixParameter(11,csigma1000);
  func->FixParameter(12,c2sigma1000);
  
  //---- 1000 multi-compton
  func->FixParameter(13,Nratio1000);
  func->FixParameter(14,gausmean1000);
  func->FixParameter(15,gaussig1000);
  
  //---- 482 bi207
  func->SetParameter(17,Ek482);
  func->SetParLimits(17,Ek482*0.8,Ek482*1.2);
  func->SetParameter(18,s482);
  func->SetParLimits(18,s482*0.5,s482*1.2);
  func->FixParameter(19,calib);
  
  //---- 976 bi207
  func->SetParameter(21,Ek976);
  func->SetParLimits(21,Ek976*0.8,Ek976*1.2);
  func->SetParameter(22,s976);
  func->SetParLimits(22,s976*0.5,s976*1.2);
  func->FixParameter(23,calib);
    
  //---- offset
  func->FixParameter(24,func->GetParameter(24));
  
  if(lowres || !thick){
    //---- 500 multi-compton
    func->FixParameter(6,0.0);
    func->FixParameter(7,1.0e-2);
    func->FixParameter(8,1.0e-2);
  }
    
  if(lowres){
    
    Float_t m=0.7;
    Float_t p=1.3;
    
    //---- 500 compton
    func->ReleaseParameter(2);
    func->SetParameter(2,N500);
    func->SetParLimits(2,0.01,1.0e4);
    func->ReleaseParameter(4);
    func->SetParameter(4,csigma500);
    func->SetParLimits(4,csigma500*m,csigma500*p);
    func->ReleaseParameter(5);
    func->SetParameter(5,c2sigma500);
    func->SetParLimits(5,c2sigma500*m,c2sigma500*p);
    
    //---- 1000 compton
    func->ReleaseParameter(9);
    func->SetParameter(9,N1000);
    func->SetParLimits(9,0.01,1.0e4);
    func->ReleaseParameter(11);
    func->SetParameter(11,csigma1000);
    func->SetParLimits(11,csigma1000*m,csigma1000*p);
    func->ReleaseParameter(12);
    func->SetParameter(12,c2sigma1000);
    func->SetParLimits(12,c2sigma1000*m,c2sigma1000*p);
    
    //---- 1000 multi-compton
    func->FixParameter(13,0.0);
    func->FixParameter(14,1.0e-2);
    func->FixParameter(15,1.0e-2);
    
    //---- 482 bi207
    //func->SetParameter(16,0.0);
    //func->SetParLimits(16,0.0,1.0e6);
      
  }
  
  ///////////////////////////////////////////////////////////////////
  TCanvas *c73 = new TCanvas("c73","- Rough Fit 2 -",xmove,0,xsize,ysize);
  c73->cd();
  xmin=280*calib;
  xmax=1200*calib;
  if(xmax>hmax) xmax=hmax;
  hbeta1->Draw();
  hbeta1->Fit("func",fitopt,"",xmin,xmax);
  ///////////////////////////////////////////////////////////////////
  
  N500         = func->GetParameter(2);
  cedge500     = func->GetParameter(3);
  csigma500    = func->GetParameter(4);
  c2sigma500   = func->GetParameter(5);
  Nratio500    = func->GetParameter(6);
  gausmean500  = func->GetParameter(7);
  gaussig500   = func->GetParameter(8);
    
  N1000        = func->GetParameter(9);
  cedge1000    = func->GetParameter(10);
  csigma1000   = func->GetParameter(11);
  c2sigma1000  = func->GetParameter(12);
  Nratio1000   = func->GetParameter(13);
  gausmean1000 = func->GetParameter(14);
  gaussig1000  = func->GetParameter(15);
    
  Ek482        = func->GetParameter(17);
  Ek482_err    = func->GetParError(17);
  s482         = func->GetParameter(18);
  s482_err     = func->GetParError(18);
  cal482       = func->GetParameter(19);
  cal482_err   = func->GetParError(19);
  
  Ek976        = func->GetParameter(21);
  Ek976_err    = func->GetParError(21);
  s976         = func->GetParameter(22);
  s976_err     = func->GetParError(22);
  cal976       = func->GetParameter(23);
  cal976_err   = func->GetParError(23);
    
  calib        = (Ek976-Ek482)/494.;
  if(lowres || thick){
    if(lowres)calib = cal976;
    if(thick) calib = (Ek976-cedge500)/(eg_diff-Eloss);
    cal482     = calib;
    cal976     = calib;
  }
  
  scint_quench = calib*((976.-Eloss)-(Ek976/calib));
  nEk976       = Ek976;
  nEk482       = Ek482;
  if(scint_quench>0){
    nEk976     = Ek976+scint_quench;
    nEk482     = Ek482+scint_quench;
  }  
  ns976        = s976;
  if(sped<s976) ns976 = TMath::Sqrt((s976*s976)-(sped*sped));
  ns482        = s482;
  if(sped<s482) ns482 = TMath::Sqrt((s482*s482)-(sped*sped));
  
  tot_chi2=func->GetChisquare();
  tot_ndf=func->GetNDF();
    
  pmin=xmin-(3*csigma500);
  pmax=xmax+(4*s976);
  if(pmin<hmin) pmin=hmin;
  if(pmax>hmax) pmax=hmax;
  hbeta1->SetAxisRange(pmin,pmax,"X");
  c73->Update();
  c73->Print("AA_bi207_rough2.png");
  ctmp->cd();
  
  cout<<"\n\n ROUGH FIT 2"<<endl;
  megaPrint();
  
  
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //      FINAL FIT TO BI207 SPECTRUM 
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
    
  //---- expo
  func->FixParameter(0,func->GetParameter(0));
  func->FixParameter(1,func->GetParameter(1));
  
  //---- 500 compton
  func->FixParameter(2,N500);
  func->FixParameter(3,cedge500);
  func->FixParameter(4,csigma500);
  func->FixParameter(5,c2sigma500);
    
  //---- 500 multi-compton
  func->FixParameter(6,Nratio500);
  func->FixParameter(7,gausmean500);
  func->FixParameter(8,gaussig500);
    
  //---- 1000 compton
  func->FixParameter(9,N1000);
  func->FixParameter(10,cedge1000);
  func->FixParameter(11,csigma1000);
  func->FixParameter(12,c2sigma1000);
    
  //---- 1000 multi-compton
  func->FixParameter(13,Nratio1000);
  func->FixParameter(14,gausmean1000);
  func->FixParameter(15,gaussig1000);
     
  //---- 482 bi207
  func->SetParameter(17,Ek482);
  func->SetParLimits(17,Ek482*0.8,Ek482*1.2);
  func->SetParameter(18,s482);
  func->SetParLimits(18,s482*0.8,s482*1.2);
  func->FixParameter(19,calib);
    
  //---- 976 bi207
  func->SetParameter(21,Ek976);
  func->SetParLimits(21,Ek976*0.8,Ek976*1.2);
  func->SetParameter(22,s976);
  func->SetParLimits(22,s976*0.8,s976*1.2);
  func->FixParameter(23,calib);
  
  //---- offset
  func->FixParameter(24,func->GetParameter(24));
  
  if(lowres || !thick){
    //---- 500 multi-compton
    func->FixParameter(6,0.0);
    func->FixParameter(7,1.0e-2);
    func->FixParameter(8,1.0e-2);
  }
  
  if(lowres){
    //---- 1000 multi-compton
    func->FixParameter(13,0.0);
    func->FixParameter(14,1.0e-2);
    func->FixParameter(15,1.0e-2);
  }
  
  ///////////////////////////////////////////////////////////////////
  TCanvas *c74 = new TCanvas("c74","- Final Fit -",xmove,0,xsize,ysize);
  c74->cd();
  xmin=280*calib;
  xmax=1200*calib;
  //xmax=Ek976+(10*s976);
  if(xmax>hmax) xmax=hmax;
  hbeta2->Draw();
  hbeta2->Fit("func",fitopt,"",xmin,xmax);
  func->DrawCopy("LSAME");
  ///////////////////////////////////////////////////////////////////
  
  N500         = func->GetParameter(2);
  cedge500     = func->GetParameter(3);
  csigma500    = func->GetParameter(4);
  c2sigma500   = func->GetParameter(5);
  Nratio500    = func->GetParameter(6);
  gausmean500  = func->GetParameter(7);
  gaussig500   = func->GetParameter(8);
    
  N1000        = func->GetParameter(9);
  cedge1000    = func->GetParameter(10);
  csigma1000   = func->GetParameter(11);
  c2sigma1000  = func->GetParameter(12);
  Nratio1000   = func->GetParameter(13);
  gausmean1000 = func->GetParameter(14);
  gaussig1000  = func->GetParameter(15);
    
  Ek482        = func->GetParameter(17);
  Ek482_err    = func->GetParError(17);
  s482         = func->GetParameter(18);
  s482_err     = func->GetParError(18);
  cal482       = func->GetParameter(19);
  cal482_err   = func->GetParError(19);
  
  Ek976        = func->GetParameter(21);
  Ek976_err    = func->GetParError(21);
  s976         = func->GetParameter(22);
  s976_err     = func->GetParError(22);
  cal976       = func->GetParameter(23);
  cal976_err   = func->GetParError(23);
  
  
  ///////////////////////////////////////////////////////////////////
  TF1 *k = new TF1("k",f_976K,Ek976-(3*s976),Ek976+(3*s976),4);
  k->FixParameter(0,func->GetParameter(20));
  k->FixParameter(1,Ek976);
  k->FixParameter(2,s976);
  k->FixParameter(3,cal976);
  k->SetLineColor(kBlue);
  k->DrawCopy("LSAME");
  TF1 *l = new TF1("l",f_976L,Ek976+(72.3*cal976)-(3*s976),Ek976+(72.3*cal976)+(3*s976),4);
  l->FixParameter(0,func->GetParameter(20));
  l->FixParameter(1,Ek976);
  l->FixParameter(2,s976);
  l->FixParameter(3,cal976);
  l->SetLineColor(kBlue+1);
  l->DrawCopy("LSAME");
  TF1 *m = new TF1("m",f_976M,Ek976+(84.3*cal976)-(3*s976),Ek976+(84.3*cal976)+(3*s976),4);
  m->FixParameter(0,func->GetParameter(20));
  m->FixParameter(1,Ek976);
  m->FixParameter(2,s976);
  m->FixParameter(3,cal976);
  m->SetLineColor(kBlue+2);
  m->DrawCopy("LSAME");
  ///////////////////////////////////////////////////////////////////
  
  
  calib        = (Ek976-Ek482)/494.;
  if(lowres || thick){
    if(lowres)calib = cal976;
    if(thick) calib = (Ek976-cedge500)/(eg_diff-Eloss);
    cal482     = calib;
    cal976     = calib;
  }
  
  scint_quench = calib*((976.-Eloss)-(Ek976/calib));
  nEk976       = Ek976;
  nEk482       = Ek482;
  if(scint_quench>0){
    nEk976     = Ek976+scint_quench;
    nEk482     = Ek482+scint_quench;
  }    
  ns976        = s976;
  if(sped<s976) ns976 = TMath::Sqrt((s976*s976)-(sped*sped));
  ns482        = s482;
  if(sped<s482) ns482 = TMath::Sqrt((s482*s482)-(sped*sped));
  
  tot_chi2=func->GetChisquare();
  tot_ndf=func->GetNDF();
    
  pmin=xmin-(3*csigma500);
  pmax=xmax+(4*s976);
  if(pmin<hmin) pmin=hmin;
  if(pmax>hmax) pmax=hmax;
  hbeta2->SetAxisRange(pmin,pmax,"X");
    
  ostringstream oss;
  oss<<betarun<<"  ==>  "<<ns976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)<<" % FWHM @ 1MeV ";
  hbeta2->SetTitle(oss.str().c_str());
  c74->Update();
  c74->Print("AA_bi207_final.png");
  c74->Print("AA_bi207_final.pdf");
  ctmp->cd();
  
  cout<<"\n\n FINAL FIT RESULTS"<<endl;
  megaPrint();
    
  delete func;
  delete ctmp;
  
  
  //---- do a fit to the subtracted histogram!
  ///////////////////////////////////////////
  quickfit_976kev(sub,Ek976-(1.3*s976),Ek976,((1060.-Eloss)*calib)+(2.0*s976),2,"AA_quickfit-976keV");
  quickfit_landau_976(sub,Ek976,s976,Eloss,2,"AA_quickfit-landau-976keV");
  
  
  return 0;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void megaPrint(){
  cout<<"======================================================================"<<endl;
  cout<<describe<<endl;
  cout<<"\t Pedestal      = "<<pedrun<<endl;
  cout<<"\t Beta File     = "<<betarun<<endl;
  cout<<"\t Gamma File    = "<<gammarun<<endl;
  
  cout<<"\n Energy Calibration from Ek976-Ek482 = "
      <<1/calib<<" keV/ch"<<endl;
  cout<<"\t Calib 976-482 = "<<calib<<endl;
  cout<<"\t Calib 482 keV = "<<cal482<<" +/- "<<cal482_err<<endl;
  cout<<"\t Calib 976 keV = "<<cal976<<" +/- "<<cal976_err<<endl;
  
  cout<<"\n 976keV peak fit "<<endl;
  cout<<"\t 976keV Mean   = "<<Ek976<<" +/- "<<Ek976_err
      <<" = "<<Ek976/cal976<<" keV"<<endl;
  cout<<"\t 976keV Sigma  = "<<s976<<" +/- "<<s976_err
      <<" = "<<s976/cal976<<" keV"<<endl;
  cout<<"\t FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "
      <<s976_err/Ek976*2.354*100<<" % "<<endl;
  cout<<"\t FWHM @ 1MeV   = "
      <<s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)
      <<" +/- "
      <<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)
      <<" % "<<endl;
  
  cout<<"\n 482keV peak fit "<<endl;
  cout<<"\t 482keV Mean   = "<<Ek482<<" +/- "<<Ek482_err
      <<" = "<<Ek482/cal482<<" keV"<<endl;
  cout<<"\t 482keV Sigma  = "<<s482<<" +/- "<<s482_err
      <<" = "<<s482/cal482<<" keV"<<endl;
  cout<<"\t FWHM @ 482keV = "<<s482/Ek482*2.354*100<<" +/- "
      <<s482_err/Ek482*2.354*100<<" % "<<endl;
  cout<<"\t FWHM @ 1MeV   = "
      <<s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)
      <<" +/- "
      <<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)
      <<" % "<<endl;
  
  cout<<"\n 976keV Scint Quench and Pedestal Sigma Contribution "<<endl;
  cout<<"\t Scint Quench  = "<<scint_quench<<" bins = "
      <<scint_quench/cal976<<" keV"<<endl;
  cout<<"\t New 976 Mean  = "<<nEk976<<" bins = "
      <<nEk976/cal976<<" keV"<<endl;
  cout<<"\t New 976 Sigma = "<<ns976<<" bins = "
      <<ns976/cal976<<" keV"<<endl;
  cout<<"\t FWHM @ 976keV = "<<ns976/nEk976*2.354*100<<" +/- "
      <<s976_err/nEk976*2.354*100<<" % "<<endl;
  cout<<"\t FWHM @ 1MeV   = "
      <<ns976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)
      <<" +/- "
      <<s976_err/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)
      <<" % "<<endl;
  
  cout<<"\n 482keV Scint Quench and Pedestal Sigma Contribution "<<endl;
  cout<<"\t Scint Quench  = "<<scint_quench<<" bins = "
      <<scint_quench/cal482<<" keV"<<endl;
  cout<<"\t New 482 Mean  = "<<nEk482<<" bins = "
      <<nEk482/cal482<<" keV"<<endl;
  cout<<"\t New 482 Sigma = "<<ns482<<" bins = "
      <<ns482/cal482<<" keV"<<endl;
  cout<<"\t FWHM @ 482keV = "<<ns482/nEk482*2.354*100<<" +/- "
      <<s482_err/nEk482*2.354*100<<" % "<<endl;
  cout<<"\t FWHM @ 1MeV   = "
      <<ns482/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)
      <<" +/- "
      <<s482_err/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)
      <<" % "<<endl;
  
  cout<<"\n FIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "
      <<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
  cout<<"======================================================================\n\n"<<endl;
  
  outfile<<"======================================================================"<<endl;
  outfile<<describe<<endl;
  outfile<<"\t Pedestal      = "<<pedrun<<endl;
  outfile<<"\t Beta File     = "<<betarun<<endl;
  outfile<<"\t Gamma File    = "<<gammarun<<endl;
  
  outfile<<"\n Energy Calibration from Ek976-Ek482 = "<<1/calib
	 <<" keV/ch"<<endl;
  outfile<<"\t Calib 976-482 = "<<calib<<endl;
  outfile<<"\t Calib 482 keV = "<<cal482<<" +/- "<<cal482_err<<endl;
  outfile<<"\t Calib 976 keV = "<<cal976<<" +/- "<<cal976_err<<endl;
  
  outfile<<"\n 976keV peak fit "<<endl;
  outfile<<"\t 976keV Mean   = "<<Ek976<<" +/- "<<Ek976_err
	 <<" = "<<Ek976/cal976<<" keV"<<endl;
  outfile<<"\t 976keV Sigma  = "<<s976<<" +/- "<<s976_err
	 <<" = "<<s976/cal976<<" keV"<<endl;
  outfile<<"\t FWHM @ 976keV = "<<s976/Ek976*2.354*100<<" +/- "
	 <<s976_err/Ek976*2.354*100<<" % "<<endl;
  outfile<<"\t FWHM @ 1MeV   = "
	 <<s976/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)
	 <<" +/- "
	 <<s976_err/Ek976*2.354*100*TMath::Sqrt(Ek976/cal976/1000)
	 <<" % "<<endl;
  
  outfile<<"\n 482keV peak fit "<<endl;
  outfile<<"\t 482keV Mean   = "<<Ek482<<" +/- "
	 <<Ek482_err<<" = "<<Ek482/cal482<<" keV"<<endl;
  outfile<<"\t 482keV Sigma  = "<<s482<<" +/- "
	 <<s482_err<<" = "<<s482/cal482<<" keV"<<endl;
  outfile<<"\t FWHM @ 482keV = "<<s482/Ek482*2.354*100
	 <<" +/- "<<s482_err/Ek482*2.354*100<<" % "<<endl;
  outfile<<"\t FWHM @ 1MeV   = "
	 <<s482/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)
	 <<" +/- "
	 <<s482_err/Ek482*2.354*100*TMath::Sqrt(Ek482/cal482/1000)
	 <<" % "<<endl;
  
  outfile<<"\n 976keV Scint Quench and Pedestal Sigma Contribution "<<endl;
  outfile<<"\t Scint Quench  = "<<scint_quench<<" bins = "
	 <<scint_quench/cal976<<" keV"<<endl;
  outfile<<"\t New 976 Mean  = "<<nEk976<<" bins = "
	 <<nEk976/cal976<<" keV"<<endl;
  outfile<<"\t New 976 Sigma = "<<ns976<<" bins = "
	 <<ns976/cal976<<" keV"<<endl;
  outfile<<"\t FWHM @ 976keV = "<<ns976/nEk976*2.354*100
	 <<" +/- "<<s976_err/nEk976*2.354*100<<" % "<<endl;
  outfile<<"\t FWHM @ 1MeV   = "
	 <<ns976/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)
	 <<" +/- "
	 <<s976_err/nEk976*2.354*100*TMath::Sqrt(nEk976/cal976/1000)
	 <<" % "<<endl;
  
  outfile<<"\n 482keV Scint Quench and Pedestal Sigma Contribution "<<endl;
  outfile<<"\t Scint Quench  = "<<scint_quench<<" bins = "
	 <<scint_quench/cal482<<" keV"<<endl;
  outfile<<"\t New 482 Mean  = "<<nEk482<<" bins = "
	 <<nEk482/cal482<<" keV"<<endl;
  outfile<<"\t New 482 Sigma = "<<ns482<<" bins = "
	 <<ns482/cal482<<" keV"<<endl;
  outfile<<"\t FWHM @ 482keV = "<<ns482/nEk482*2.354*100
	 <<" +/- "<<s482_err/nEk482*2.354*100<<" % "<<endl;
  outfile<<"\t FWHM @ 1MeV   = "
	 <<ns482/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)
	 <<" +/- "
	 <<s482_err/nEk482*2.354*100*TMath::Sqrt(nEk482/cal482/1000)
	 <<" % "<<endl;
  
  outfile<<"\n FIT RESULTS, chi2/ndf = "<<tot_chi2<<" / "
	 <<tot_ndf<<" = "<<tot_chi2/tot_ndf<<endl;
  outfile<<"======================================================================\n\n"<<endl;
  
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

