/***********************************************************
A MACRO TO FIT THE 976KEV PEAK WITH 3 GAUSSIANS
------------------------------------------------------------
This program fits the 976keV Bi207 peak using 3 gaussians
for the KLM shells.

VERSION:   09.09.23

change log: (most recent at top)
-----------------------------------
23sep ~ TH1F *histio --> TH1 *histio
      + resolution at 1MeV
19jun ~ tidy up a bit
      + added the scint quenching factor of 30keV
        when calculating the calibration

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/
#include "./include/fittingFuncs.hxx"
#include "./include/manip_histos.hxx"

void quickfit_976kev(TH1 *histio, Int_t fmin=880, Int_t fmid=960, Int_t fmax=1140, Int_t rebin=1, string name=""){
  cout<<"\n-------------------------------------------------------------------------------------\n";
  cout<<  "USAGE ==> quickfit_976kev.cxx (hist_ptr, fit_min, fit_mid, fit_max, rebin, plot name)\n";
  cout<<  "-------------------------------------------------------------------------------------\n\n";
  
  beautify();
  
  Float_t Eloss=17.0;  // keV
  Float_t quench=30.0; //in keV
  //quench = 0.0;  // set to 0 for simulations?
  
  TCanvas *bcan = new TCanvas("bcan","976 klm fit",0,0,700,500);
  bcan->cd();
  TH1F *bhist =(TH1F*)cloneHist((TH1F*)histio,"bhist","Fit 976keV with KLM Gaussians");
  bhist->Rebin(rebin);
  
  TF1 *pfind=new TF1("pfind","gaus",fmin,fmax);
  bhist->Fit("pfind","R");
  Float_t mean=pfind->GetParameter(1);
  Float_t sig=pfind->GetParameter(2);
  bhist->Fit("pfind","","",mean-sig,mean+sig);
  mean=pfind->GetParameter(1);
  //Float_t calib=mean/(976.-Eloss);
  Float_t calib=mean/(976.-quench-Eloss);
  
  TF1 *bfunc = new TF1("bfunc",f_976kev,fmin,fmax,4);
  bfunc->SetLineColor(2);
  bfunc->SetLineWidth(2);
  //bfunc->SetParNames("Normal","EnergyK","SigmaK","Calib");
  bfunc->SetParNames("Gn976","Ge976","Gs976","Gc976");
  
  bfunc->SetParameter(0,1e3);     // 976 keV normal
  bfunc->SetParLimits(0,1,1e8);
  bfunc->SetParameter(1,mean);    // 976 keV peak
  bfunc->SetParLimits(1,1,3000);
  bfunc->SetParameter(2,sig);     // 976 keV sigma
  bfunc->SetParLimits(2,1,500);
  bfunc->SetParameter(3,calib);   // 976 calib
  bfunc->SetParLimits(3,0.1,2.5);
  
  // SHOULD WE FIX THESE?
  bfunc->FixParameter(3,calib);  // FIX THE CALIBRATION?
  
  bhist->Fit("bfunc","R");
    
  Float_t resol=(bfunc->GetParameter(2)/bfunc->GetParameter(1))*2.35*100.;
  cout<<"\n\n\t Resolution    =  "<<resol<<" % \n";
  Float_t resmev=resol*TMath::Sqrt(bfunc->GetParameter(1)/bfunc->GetParameter(3)/1000);
  cout<<"    \t Res @ 1MeV    =  "<<resmev<<" % \n";
  Float_t energy=bfunc->GetParameter(1)/bfunc->GetParameter(3);
  cout<<    "\t 976keV Energy =  "<<energy<<" keV \n";
  cout<<    "\t Energy Loss   =  "<<Eloss<<" keV \n";
  cout<<"\n\n";
  
  Int_t resol2=Int_t(resol*100);
  Int_t resmev2=Int_t(resmev*100);
  Int_t energ2=Int_t(energy*10);
  ostringstream oss_resol,oss_energy,oss_eloss;
  //oss_resol<<"Resolution  =  "<<resol2/100.<<" % ";
  oss_resol<<"Res @ 1MeV  =  "<<resmev2/100.<<" % ";
  oss_energy<<"Energy  =  "<<energ2/10.<<" keV ";
  oss_eloss<<"E-loss  =  "<<Eloss<<" keV ";
  
  bhist->SetAxisRange(fmin*0.8,fmax*1.4,"X");
  Double_t xpos=bfunc->GetParameter(1)+(4*bfunc->GetParameter(2));
  Double_t ypos=bhist->GetBinContent(bhist->GetMaximumBin());
  Float_t fsize=0.045; Color_t color=kBlue;
  TLatex *tex_v1 = new TLatex(xpos,ypos*0.6,oss_resol.str().c_str());
  tex_v1->SetTextSize(fsize); tex_v1->SetTextColor(color);
  tex_v1->Draw();
  TLatex *tex_v2 = new TLatex(xpos,ypos*0.5,oss_energy.str().c_str());
  tex_v2->SetTextSize(fsize); tex_v2->SetTextColor(color);
  tex_v2->Draw();
  TLatex *tex_v3 = new TLatex(xpos,ypos*0.4,oss_eloss.str().c_str());
  tex_v3->SetTextSize(fsize); tex_v3->SetTextColor(color);
  tex_v3->Draw();
  
  bcan->Update();
  if(name==""){
    bcan->Print("quickfit-976kev.png");
  }else{  
    bcan->Print((name+".png").c_str());
    bcan->Print((name+".pdf").c_str());
  }
  
  
  delete pfind;
  delete bfunc;
}
