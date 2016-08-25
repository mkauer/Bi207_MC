gROOT->Reset();

#include "./include/manip_files.hxx"
#include "./include/manip_histos.hxx"

///////////////////////////////////////////////////////////
const Int_t input_type = 1;  // 0==rootfile : 1==datafile
const Int_t bins       = 4000;
const Int_t histmax    = 4000;
///////////////////////////////////////////////////////////

Int_t _combineMC(string fname="combined", string outdir="./DBtemp/")
{
  beautify();
  
  if(input_type==0){
    char betafile[256]="MCbetas/results.root";
    char gammafile[256]="MCgammas/results.root";
    char hname[24]="h100";
    
    TFile *file1=new TFile(betafile,"READ");
    TH1F *MCbeta=(TH1F*)file1->Get(hname);
    MCbeta->SetName("MCbeta");
    MCbeta->SetTitle("MCbeta");
    
    TFile *file2=new TFile(gammafile,"READ");
    TH1F *MCgamma=(TH1F*)file2->Get(hname);
    MCgamma->SetName("MCgamma");
    MCgamma->SetTitle("MCgamma");
  }
  
  if(input_type==1){
    char betafile[256]="MCbetas/output.dat";
    char gammafile[256]="MCgammas/output.dat";
    
    TH1F *MCbeta = readfile(betafile,0,0,2,"MCbeta",1,histmax);
    TH1F *MCgamma = readfile(gammafile,0,0,2,"MCgamma",1,histmax);
    MCbeta->SetTitle("MCbeta");
    MCgamma->SetTitle("MCgamma");
  }
  
  TCanvas *look=new TCanvas("look","look",0,0,800,600);
  look->Divide(1,2,.02,.02,0);
  look->cd(1);
  MCbeta->Draw();
  look->Update();
  look->cd(2);
  MCgamma->Draw();
  look->Update();
  
  TString rootname=outdir+fname;
  string s1=rootname+".png";
  string s2=rootname+".C";
  look->Print(s1.c_str());
  look->Print(s2.c_str());
  
  TFile *outfile=gROOT->FindObject(rootname+".root");
  if(outfile) outfile->Close();
  outfile=new TFile(rootname+".root","RECREATE");
  MCbeta->Write();
  MCgamma->Write();
  outfile->Close();
  
  if(genConf(fname,outdir)) return 1;
  
  cout<<"\n\n\tDONE CREATING ==> "<<rootname+".root"<<" \n\n\n";

  return 0;
}


//--- GENERATE A GENERAL CONFIG FILE FOR THE RUNS
Int_t genConf(string fname="combined", string outdir="./DBtemp/"){
  ofstream fout;
  fout.open((outdir+fname+".conf").c_str());
  
  fout<<"void config(){ \n\n";
  
  fout<<"// example = 8in-Ham-SBA-1900v__25x22x10-hex-ej200-esr-mylar__bi-center-6cm \n";
  fout<<"// example = 8in-ETL-1500v_LG-pol-tef__14x14x2-bc404-esr__bi-center-6cm \n\n";
  
  //fout<<"     describe  = \"pmt-info__scint-info__source-info\"; \n\n";
  fout<<"     describe  = \""<<fname.c_str()<<"\"; \n\n";
  
  fout<<"     pedrun    = \"\";          // pedestal run-name \n";
  fout<<"     betarun   = \"MCbeta\";    // beta+gamma run-name \n";
  fout<<"     gammarun  = \"MCgamma\";   // gamma run-name \n\n";
  
  fout<<"     gamma     = 400;     // the 400 kev COMPTON edge position \n";
  fout<<"     beta      = 976;     // the 976 kev CONVERSION-E peak position \n";
  fout<<"     Eloss     = 17.0;    // energy loss in air and reflector \n\n";
  
  fout<<"     lowres    = 0;       // if resolution is >8% set this true \n";
  fout<<"     thick     = 1;       // if using thick scint, set this true \n";
  fout<<"     fitopt    = \"\";      // global fit options (chi^2 or likelyhood etc) \n";
  fout<<"     rebin     = 2;       // rebinning of the histrograms \n";
  fout<<"     maxbin    = "<<histmax<<";    // maximum range you want in your histos \n\n";
  
  fout<<"} \n";
  
  fout.close();  

  return 0;
}

