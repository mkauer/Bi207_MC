#include "../include/manip_files.hxx"
#include "../include/manip_histos.hxx"

Int_t _histoMC(char input[64]="output.dat",Int_t bins=4000,Int_t max=4000)
{
  beautify();
    
  TH1F *histo = readfile(input,0,0,2,"histo",1,max);
  
  TCanvas *c10 = new TCanvas("c10","c10",700,500);
  histo->Draw();
  c10->Update();
  c10->Print("_histo.png");
  
  return 0;
}

