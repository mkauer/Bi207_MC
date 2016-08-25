/***********************************************************
A MACRO TO FIT THE 976KEV PEAK WITH LANDAU ENERGY LOSSES
------------------------------------------------------------
This program fits the 976keV Bi207 peak using 3 gaussians
for the KLM shells and also takes into account landau losses
for the mean energy loss you specify.

VERSION: 09.09.23

change log: (most recent at top)
-----------------------------------
23sep ~ TH1F *histio --> TH1 *histio
      + resolution at 1MeV
19jul ~ tried fixing the calibration and engery
        loss parameters and seems to work well.
      + added the scint quenching factor of 30keV
        when calculating the calibration

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/
#include "./include/fittingFuncs.hxx"
#include "./include/manip_histos.hxx"

void quickfit_landau_976(TH1 *histio, Double_t Ge976=960, Double_t Gs976=30, Double_t Eloss=17, Int_t rebin=1, string name=""){
  cout<<"-----------------------------------------------------------------------------------\n";
  cout<<"USAGE ==> quickfit_landau_976.cxx (hist_ptr, Ge976, Gs976, Eloss, rebin, plot name)\n";
  cout<<"-----------------------------------------------------------------------------------\n";
  
  beautify();
  
  Float_t quench=30.0; //in keV
  //quench = 0.0;  // set to 0 for simulations?
  
  TCanvas *lcan = new TCanvas("lcan","976 landau fit",0,0,700,500);
  TH1F *lhist =(TH1F*)cloneHist((TH1F*)histio,"lhist","Fit 976keV with Landau Tail");
  lhist->Rebin(rebin);
  
  Int_t fmin=Ge976-(5*Gs976);
  //Int_t fmax=Ge976+(6*Gs976);
  Int_t fmax=Ge976+(3*Gs976);
  Int_t pmin=Ge976-(10*Gs976);
  Int_t pmax=Ge976+(18*Gs976);
  
  Float_t calib=Ge976/(976.-quench-Eloss);
  
  TF1 *lfunc = new TF1("lfunc",f_landau_976,fmin,fmax,8);
  lfunc->SetLineColor(2);
  lfunc->SetLineWidth(2);
  lfunc->SetParNames("Gn976","Ge976","Gs976","Gc976","Ln976","Eloss","Ls976","offset");
  
  lfunc->SetParameter(0,1e3);    // 976 keV normal
  lfunc->SetParLimits(0,0,1e7);
  lfunc->SetParameter(1,Ge976);  // 976 keV peak
  lfunc->SetParLimits(1,Ge976*.5,Ge976*1.5);
  lfunc->SetParameter(2,Gs976);  // 976 keV sigma
  lfunc->SetParLimits(2,Gs976*.5,Gs976*1.5);
  lfunc->SetParameter(3,calib);  // 976 calib
  lfunc->SetParLimits(3,.5,3);
  lfunc->SetParameter(4,1e3);    // 976 landau normal
  lfunc->SetParLimits(4,0,1e7);
  lfunc->SetParameter(5,Eloss);  // 976 mean energy losses
  lfunc->SetParLimits(5,1,50);
  lfunc->SetParameter(6,Gs976);  // 976 landau sigma
  lfunc->SetParLimits(6,0,500);
  lfunc->SetParameter(7,1);      // the offset
  lfunc->SetParLimits(7,0,500);
  
  // testing
  //lfunc->SetParameter(8,1.851);  // L-shell normalization
  //lfunc->SetParLimits(8,1.840,2.001);
  //lfunc->SetParameter(9,72.3);   // L-shell shift from K-shell
  //lfunc->SetParLimits(9,60.1,85.1);
  
  // SHOULD WE FIX THESE?
  lfunc->FixParameter(3,calib);  // FIX THE CALIBRATION?
  lfunc->FixParameter(5,Eloss);  // FIX THE ENERGY LOSSES?
  
  lhist->Fit("lfunc","","",fmin,fmax);
  
  Float_t resol=(lfunc->GetParameter(2)/lfunc->GetParameter(1))*2.35*100;
  cout<<"\n\n\t Resolution    =  "<<resol<<" % \n";
  Float_t resmev=resol*TMath::Sqrt(lfunc->GetParameter(1)/lfunc->GetParameter(3)/1000);
  cout<<"    \t Res @ 1MeV    =  "<<resmev<<" % \n";
  Float_t energy=lfunc->GetParameter(1)/lfunc->GetParameter(3);
  cout<<    "\t 976keV Energy =  "<<energy<<" keV \n";
  Eloss=lfunc->GetParameter(5);
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
  
  lhist->SetAxisRange(fmin*0.8,fmax*1.4,"X");
  Double_t xpos=lfunc->GetParameter(1)+(3*lfunc->GetParameter(2));
  Double_t ypos=lhist->GetBinContent(lhist->GetMaximumBin());
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
  
  lcan->Update();
  if(name==""){
    lcan->Print("quickfit-landau-976kev.png");
  }else{
    lcan->Print((name+".png").c_str());
    lcan->Print((name+".pdf").c_str());
  }
  
  
  delete lfunc;
}

