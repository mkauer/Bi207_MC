/***********************************************************
FUNCTIONS FOR MANIPULATING HISTOGRAMS
------------------------------------------------------------

VERSION: 09.08.31

change log: (most recent at top)
-----------------------------------
31aug + include init_memory
19jul ~ changed beautify a little
16jan ~ added shift==0 to shiftHist
16jan ~ typo changed in drawMix()
04nov + function to return pointer to subtracted hist
04nov + function to draw subtracted histo
02nov + function to normalize two histos to each other
31oct + function for shifting a histo (+) or (-) in x-axis
29oct + function to resize a histogram

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/

#include "init_memory.cxx"

//--- SET ROOT DEFAULTS TO LOOK GOOD
void beautify()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadLeftMargin  (0.10);
  gStyle->SetPadTopMargin   (0.08);
  gStyle->SetPadRightMargin (0.15);
  
  gStyle->SetOptStat("RMen");
  gStyle->SetOptFit(112);
  
  gStyle->SetFillColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetStatW(0.13);
  gStyle->SetStatH(0.08);
}

//--- RETURN A CLONED HISTOGRAM
TH1F* cloneHist(TH1F *inhist,char *name="Name",char *title="Title")
{
  Int_t nbins=inhist->GetXaxis()->GetNbins();
  Int_t hmin=inhist->GetXaxis()->GetBinLowEdge(1);
  Int_t hmax=inhist->GetXaxis()->GetBinUpEdge(nbins);
  Double_t content;
  Double_t nevents=inhist->GetEntries();
  
  //cout<<"\n\n\t nbins = "<<nbins<<"  hmin = "<<hmin<<"  hmax = "<<hmax<<"\n\n";
  
  TH1F *outhist=new TH1F(name,title,nbins,hmin,hmax);
  
  for(Int_t i=0;i<nbins;i++){
    content=inhist->GetBinContent(i);
    outhist->SetBinContent(i,content);
  }  
  outhist->SetEntries(nevents);

  return outhist;
}

//--- SHIFT A HISTOGRAM
Int_t shiftHist(TH1F *hist, Double_t double_shift)
{
  Int_t bins=hist->GetXaxis()->GetNbins();
  Double_t xmin=hist->GetXaxis()->GetBinLowEdge(1);
  Double_t xmax=hist->GetXaxis()->GetBinUpEdge(bins);
  double_shift=double_shift*(bins/(xmax-xmin));
  Int_t shift;
  if(double_shift<0) shift=TMath::FloorNint(double_shift);
  if(double_shift>0) shift=TMath::CeilNint(double_shift);
  if(shift==0) return 0;
  if(shift>0){
    for(Int_t i=1; i<=bins; i++){
      if(i+shift<=bins) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift>bins) hist->SetBinContent(i,0);
    }
    return 0;
  }
  if(shift<0){
    for(Int_t i=bins; i>0; i--){
      if(i+shift>0) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift<=0) hist->SetBinContent(i,0);
    }    
    return 0;
  }
  return 1;
} 

//--- RETURN A POINTER TO A RESIZED HISTOGRAM 
TH1F* resizeHist(TH1F *inhist,Int_t nhmin,Int_t nhmax,char *name="Name",char *title="Title")
{
  Int_t bins=inhist->GetXaxis()->GetNbins();
  Float_t hmin=inhist->GetXaxis()->GetBinLowEdge(1);
  Float_t hmax=inhist->GetXaxis()->GetBinUpEdge(bins);
  Int_t content;
  Int_t nevents=0;
  
  cout<<"\n\t bins = "<<bins<<"  hmin = "<<hmin<<"  hmax = "<<hmax<<"\n";
  
  Float_t unit=hmax/bins;
  Int_t lobin,hibin,nbins;
  if(nhmin!=0) lobin=TMath::CeilNint(nhmin/unit);
  hibin=TMath::FloorNint(nhmax/unit);
  nbins=hibin-lobin;
  
  cout<<"\t nbins = "<<nbins<<"  nhmin = "<<nhmin<<"  nhmax = "<<nhmax<<"\n";
  
  TH1F *outhist=new TH1F(name,title,nbins,nhmin,nhmax);
  
  
  for(Int_t i=lobin+1;i<=hibin;i++){
    content=inhist->GetBinContent(i-lobin);
    outhist->SetBinContent(i-lobin,content);
    nevents+=content;
  }  
  outhist->SetBinContent(0,0);
  outhist->SetBinContent(hibin+1,0);
  outhist->SetEntries(nevents);
  
  return outhist;
}

//--- NORMALIZE TWO HISTOGRAMS
void normalizer(TH1F *hgamma, TH1F *hbeta, Int_t normv=0)
{
  cout<<"\n==========================================================="<<endl;
  cout<<"Normalizing the Gamma and Beta histograms..."<<endl;
  
  Float_t beta_area,gamma_area,area_ratio;
  
  //Float_t imin=620*calib;
  //Float_t imax=800*calib;
  Float_t imin=600*calib;
  Float_t imax=700*calib;
  
  Int_t gamma_binmin=hgamma->FindBin(imin);
  Int_t gamma_binmax=hgamma->FindBin(imax);
  cout<<"\n\t Integrate Gamma between  ==>  "<<imin<<" - "<<imax<<endl;
  cout<<  "\t equals bin numbers       ==>  "<<gamma_binmin<<" - "<<gamma_binmax<<endl;
  if(normv==0) gamma_area = hgamma->TH1::Integral(gamma_binmin,gamma_binmax,"width");
  if(normv==1) gamma_area = hgamma->Integral(gamma_binmin,gamma_binmax);
  cout<<  "\t Area under gamma hist    ==>  "<<gamma_area<<endl;
  
  Int_t beta_binmin=hbeta->FindBin(imin);
  Int_t beta_binmax=hbeta->FindBin(imax);
  cout<<"\n\t Integrate Beta between   ==>  "<<imin<<" - "<<imax<<endl;
  cout<<  "\t equals bin numbers       ==>  "<<beta_binmin<<" - "<<beta_binmax<<endl;
  if(normv==0) beta_area = hbeta->TH1::Integral(beta_binmin,beta_binmax,"width");
  if(normv==1) beta_area = hbeta->Integral(beta_binmin,beta_binmax);
  cout<<  "\t Area under beta hist     ==>  "<<beta_area<<endl;
  
  area_ratio=1;
  if(Int_t(beta_area)!=0 && Int_t(gamma_area)!=0) area_ratio = beta_area/gamma_area;
  cout<<"\n\t Gamma scale factor       ==>  "<<area_ratio<<" \n"<<endl;
  if(normv==0) hgamma->Scale(area_ratio);
  if(normv==1) hgamma->TH1::Scale(area_ratio,"width");
  
  cout<<"===========================================================\n"<<endl;
}

// DRAW A CANVAS DISPLAYING THE OVERLAP OF TWO HISTOS
Int_t drawMix(TH1F *hgamma, TH1F *hbeta)
{
  TCanvas *c_mix = new TCanvas("c_mix","- Mixed -",xsize,0,xsize,ysize);
  c_mix->cd();
  TH1F *mix1 = (TH1F*)cloneHist(hgamma,"mix","Betas=blue & Gammas=green");
  mix1->SetLineColor(8);
  TH1F *mix2 = (TH1F*)cloneHist(hbeta,"mix","Betas=blue & Gammas=green");
  mix2->SetLineColor(9);
  mix1->SetAxisRange(200*calib,1400*calib,"X");
  mix2->SetAxisRange(200*calib,1400*calib,"X");
  if(mix2->GetBinContent(mix2->GetMaximumBin()) > mix1->GetBinContent(mix1->GetMaximumBin())){
    mix2->Draw();
    mix1->Draw("same");
    c_mix->Update();
    c_mix->Print("AA_mixed.png");
    return 0;
  }
  if(mix1->GetBinContent(mix1->GetMaximumBin()) > mix2->GetBinContent(mix2->GetMaximumBin())){
    mix1->Draw();
    mix2->Draw("same");
    c_mix->Update();
    c_mix->Print("AA_mixed.png");
    return 0;
  }
  return 1;
}

// DRAW A CANVAS DISPLAYING THE SUBTRACTION OF TWO HISTOS
void drawSub(TH1F *hgamma, TH1F *hbeta, Int_t rebin=1)
{
  TCanvas *c_sub = new TCanvas("c_sub","- Subtraction -",0,0,xsize,ysize);
  c_sub->cd();
  TH1F *sub1 = (TH1F*)cloneHist(hgamma,"sub","Beta - Gamma Spectra");
  TH1F *sub2 = (TH1F*)cloneHist(hbeta,"sub","Beta - Gamma Spectra");
  sub1->Rebin(rebin);
  sub2->Rebin(rebin);
  sub2->Add(sub1,-1);
  sub2->Draw();
  sub2->SetAxisRange(200*calib,1400*calib,"X");
  c_sub->Update();
  c_sub->Print("AA_subtracted.C");
  c_sub->Print("AA_subtracted.png");
}

// RETURN A POINTER TO THE SUBTRACTED HISTO
TH1F* subtract(TH1F *hgamma, TH1F *hbeta, Int_t rebin=1)
{
  TH1F *sub1 = (TH1F*)cloneHist(hgamma,"sub","Beta - Gamma Spectra");
  TH1F *sub2 = (TH1F*)cloneHist(hbeta,"sub","Beta - Gamma Spectra");
  sub1->Rebin(rebin);
  sub2->Rebin(rebin);
  sub2->Add(sub1,-1);
  return sub2;
}

