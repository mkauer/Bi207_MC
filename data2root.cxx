/***********************************************************
A FUNCTION TO READ IN RAW DATA INTO A ROOTFILE
------------------------------------------------------------
should work for many arbitrary formats of datafiles, ie the
number of columns of data, or the number of bins needed to
created the histogram.

VERSION: 09.07.19

change log: (most recent at top)
-----------------------------------
19jun ~ tidy up a bit
04jun + config file automatically generated now but you
        still must enter in the 400keV compton edge and the
        976keV CE peak in ADC bins

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/
gROOT->Reset();

#include "./include/manip_histos.hxx"
#include "./include/manip_files.hxx"


// USER VARIABLES
/////////////////////////////////////////////////////////////
Int_t numchans  = 16;   // 16 for vme and 12 for camac
Float_t histmax = 4000; // 4000 for vme and 2000 for camac
/////////////////////////////////////////////////////////////


// These parameters might need to be changed if you are trying to
// fit simulated data and the binning structure is not the same
// as the VME or CAMAC binning structure.
Float_t histmin = 0;
Int_t   binmax  = 4000; // should be the same as histmax by default
Float_t datmax  = 4000; // should be the same as histmax by default

void data2root(string file2, string file3, Int_t chan=1, string workdir="./", string outdir="./DBtemp/")
{
  beautify();
  cout<<"\n====================================================================="<<endl;
  cout<<"USAGE: vme_to_root(\"beta-file\",\"gamma-file\",channel)"<<endl;
  cout<<"USAGE: vme_to_root(\"ped-file\",\"beta-file\",\"gamma-file\",channel)"<<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  
  string file1="";
  
  ifstream infile2;
  infile2.open((workdir+file2).c_str());
  ifstream infile3;
  infile3.open((workdir+file3).c_str());
  
  if(!infile2.is_open()){
    cout<<"********** COULD NOT FIND THE FILE ==> "<<file2<<endl<<endl;
    break;
  }
  if(!infile3.is_open()){
    cout<<"********** COULD NOT FIND THE FILE ==> "<<file3<<endl<<endl;
    break;
  }
  
  cout<<"FOUND ALL THE FILES SUCCESSFULLY!"<<endl;
  
  Int_t cols2=getCol((workdir+file2).c_str());
  Int_t cols3=getCol((workdir+file3).c_str());
  
  cout<<"BETA FILE HAS THIS MANY COLUNMS  ==> "<<cols2<<endl;
  cout<<"GAMMA FILE HAS THIS MANY COLUNMS ==> "<<cols3<<endl;
  
  if(cols2!=cols3){
    cout<<"********** FILES DID NOT HAVE SAME NUMBER OF COLUNMS **********"<<endl;
    break;
  }
  
  TString rootoutfile=outdir+file2;
  TFile *outfile=gROOT->FindObject(rootoutfile+".root");
  if(outfile) outfile->Close();
  outfile=new TFile(rootoutfile+".root","RECREATE");
  
  TH1F *beta=new TH1F(file2.c_str(),file2.c_str(),binmax,histmin,histmax);
  TH1F *gamma=new TH1F(file3.c_str(),file3.c_str(),binmax,histmin,histmax);
  
  Int_t j;
  Float_t n;	      
  
  cout<<"READING IN FILE ==> "<<file2<<" ==> beta "<<endl;
  j=0;
  n=0;	      
  while(infile2>>n){
    if(j==numchans) j=0; 
    if(j==chan) beta->Fill((histmax/datmax)*n);
    //if(j==numchans) j=0; 
    j++;		
  }
  
  TCanvas *look=new TCanvas("look","look",0,0,800,600);
  look->Divide(1,2,.02,.02,0);
  look->cd(1);
  beta->Draw();
  look->Update();
  
  cout<<"READING IN FILE ==> "<<file3<<" ==> gamma "<<endl;
  j=0;
  n=0;	      
  while(infile3>>n){
    if(j==numchans) j=0; 
    if(j==chan) gamma->Fill((histmax/datmax)*n); 
    //if(j==numchans) j=0; 
    j++;		
  }
  
  look->cd(2);
  gamma->Draw();
  look->Update();
  string canvas1=outdir+file2+".png";
  string canvas2=outdir+file2+".C";
  look->Print(canvas1.c_str());
  look->Print(canvas2.c_str());
  
  beta->Write();
  gamma->Write();
 
  infile2.close();
  infile3.close();
  
  gSystem->Sleep(2000);
  gROOT->Delete("look");
  delete look;
  outfile->Close("R");
  
  genConf(file1,file2,file3,outdir);
  
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"CREATED A ROOTFILE CALLED ==>  "<<TString(rootoutfile+".root")<<endl;
  cout<<"=====================================================================\n"<<endl;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void data2root(string file1, string file2, string file3, Int_t chan=1, string workdir="./", string outdir="./DBtemp/")
{
  beautify();
  cout<<"\n====================================================================="<<endl;
  cout<<"USAGE: vme_to_root(\"beta-file\",\"gamma-file\",channel)"<<endl;
  cout<<"USAGE: vme_to_root(\"ped-file\",\"beta-file\",\"gamma-file\",channel)"<<endl;
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  
  ifstream infile1;
  infile1.open((workdir+file1).c_str());
  ifstream infile2;
  infile2.open((workdir+file2).c_str());
  ifstream infile3;
  infile3.open((workdir+file3).c_str());
  
  if(!infile1.is_open()){
    cout<<"********** COULD NOT FIND THE FILE ==> "<<file1<<endl<<endl;
    break;
  }
  if(!infile2.is_open()){
    cout<<"********** COULD NOT FIND THE FILE ==> "<<file2<<endl<<endl;
    break;
  }
  if(!infile3.is_open()){
    cout<<"********** COULD NOT FIND THE FILE ==> "<<file3<<endl<<endl;
    break;
  }
  
  cout<<"FOUND ALL THE FILES SUCCESSFULLY!"<<endl;
  
  Int_t cols1=getCol((workdir+file1).c_str());
  Int_t cols2=getCol((workdir+file2).c_str());
  Int_t cols3=getCol((workdir+file3).c_str());
  
  cout<<"PED FILE 1 HAS THIS MANY COLUNMS   ==> "<<cols1<<endl;
  cout<<"BETA FILE 2 HAS THIS MANY COLUNMS  ==> "<<cols2<<endl;
  cout<<"GAMMA FILE 3 HAS THIS MANY COLUNMS ==> "<<cols3<<endl;
  
  if(cols1!=cols2 || cols2!=cols3){
    cout<<"********** FILES DID NOT HAVE SAME NUMBER OF COLUNMS **********"<<endl;
    break;
  }
  
  TString rootoutfile=outdir+file2;
  TFile *outfile=gROOT->FindObject(rootoutfile+".root");
  if(outfile) outfile->Close();
  outfile=new TFile(rootoutfile+".root","RECREATE");
  
  Int_t binmax=4000;
  
  TH1F *ped   = new TH1F(file1.c_str(),file1.c_str(),binmax,histmin,histmax);
  TH1F *beta  = new TH1F(file2.c_str(),file2.c_str(),binmax,histmin,histmax);
  TH1F *gamma = new TH1F(file3.c_str(),file3.c_str(),binmax,histmin,histmax);
  
  Int_t numchans=cols2;
  Int_t j;
  Float_t n;	      
  
  cout<<"READING IN FILE ==> "<<file1<<" ==> ped "<<endl;
  j=0;
  n=0;
  while(infile1>>n){
    if(j==numchans)j=0; 
    if(j==chan) ped->Fill((histmax/datmax)*n,1);  
    j++;		
  }
  
  TCanvas *look=new TCanvas("look","look",0,0,1000,600);
  look->Divide(2,2,0.02,0.02,0);
  look->cd(1);
  ped->Draw();
  look->Update();
  
  cout<<"READING IN FILE ==> "<<file2<<" ==> beta "<<endl;
  j=0;
  n=0;	      
  while(infile2>>n){
    if(j==numchans) j=0; 
    if(j==chan) beta->Fill((histmax/datmax)*n,1);  
    j++;		
  }
  
  look->cd(2);
  beta->Draw();
  look->Update();
  
cout<<"READING IN FILE ==> "<<file3<<" ==> gamma "<<endl;
  j=0;
  n=0;	      
  while(infile3>>n){
    if(j==numchans) j=0; 
    if(j==chan) gamma->Fill((histmax/datmax)*n,1);  
    j++;		
  }
  
  look->cd(4);
  gamma->Draw();
  look->Update();
  string canvas1=outdir+file2+".png";
  string canvas2=outdir+file2+".C";
  look->Print(canvas1.c_str());
  look->Print(canvas2.c_str());
  
  ped->Write();
  beta->Write();
  gamma->Write();
 
  infile1.close();
  infile2.close();
  infile3.close();
  
  gSystem->Sleep(2000);
  gROOT->Delete("look");
  delete look;
  outfile->Close("R");
  
  genConf(file1,file2,file3,outdir);
  
  cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
  cout<<"CREATED A ROOTFILE CALLED ==>  "<<TString(rootoutfile+".root")<<endl;
  cout<<"=====================================================================\n"<<endl;
}


//--- GENERATE A GENERAL CONFIG FILE FOR THE RUNS
void genConf(string file1, string file2, string file3, string outdir="./DBtemp/"){
  ofstream fout;
  fout.open((outdir+file2+".conf").c_str());
  
  fout<<"void config(){ \n\n";
  
  fout<<"// example = 8in-Ham-SBA-1900v__25x22x10-hex-ej200-esr-mylar__bi-center-6cm \n";
  fout<<"// example = 8in-ETL-1500v_LG-pol-tef__14x14x2-bc404-esr__bi-center-6cm \n\n";
  
  fout<<"     describe  = \"pmt-info__scint-info__source-info\"; \n\n";
  
  fout<<"     pedrun    = \""<<file1<<"\";   // pedestal run-name \n";
  fout<<"     betarun   = \""<<file2<<"\";   // beta+gamma run-name \n";
  fout<<"     gammarun  = \""<<file3<<"\";   // gamma run-name \n\n";
  
  fout<<"     gamma     = 400;     // the 400 kev COMPTON edge position \n";
  fout<<"     beta      = 976;     // the 976 kev CONVERSION-E peak position \n";
  fout<<"     Eloss     = 14.0;    // energy loss in air and reflector \n\n";
  
  fout<<"     lowres    = 0;       // if resolution is >8% set this true \n";
  fout<<"     thick     = 0;       // if using thick scint, set this true \n";
  fout<<"     fitopt    = \"\";      // global fit options (chi^2 or likelyhood etc) \n";
  fout<<"     rebin     = 1;       // rebinning of the histrograms \n";
  fout<<"     maxbin    = "<<histmax<<";    // maximum range you want in your histos \n\n";
  
  fout<<"} \n";
  
  fout.close();  
}

