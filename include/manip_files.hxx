/***********************************************************
FUNCTIONS FOR MANIPULATING GENERAL DATA FILES
------------------------------------------------------------

VERSION: 09.07.24

change log: (most recent at top)
-----------------------------------
24jul + added the user defined part
      + added error is chan > numchans
04nov ~ this will become obsolete soon as i'm only using
        rootfiles now!
02nov ~ readfile functin already here, just rebuilding the class

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/

//--- READIN SOME DATAS FROM EITHER VME OR CAMAC DAT FILE
TH1F* readfile(char *source, Float_t ped=0, Int_t chan=1, Int_t daq=0, char *name="Name", Int_t tnumchans=16, Int_t tbinmax=4000)
{
  // VME Part
  if(daq==0){
    const Int_t numchans=16;  
    const Int_t binmax=4000;
    if(chan>=numchans){
      cout<<"\nERROR: requested channel is higher than total channels available! \n\n";
      break;
    }
    TH1F *hadc=new TH1F(name,source,binmax,0,binmax);
    ifstream infile;
    infile.open(source);
    if(infile.is_open()){
      infile.clear();
      infile.seekg(0,ios::beg);
      cout<<"\n\tREADING IN FILE ---> "<<source<<"\n";
      Int_t adc[numchans];
      Int_t j=0, events;  
      Float_t n=0;	      
      while(infile>>n){
	if(j==numchans) j=0; 
	if(j==chan) hadc->Fill(n-ped); 
	j++;		
      }
      cout<<"\n\tDONE READING FILE ---> "<<source<<"\n\n";
      infile.close();  
    }
    else cout<<"\n********** Could Not Read in File --> "<<source<<endl;
  }
  
  // CAMAC Part
  if(daq==1){ 
    const Int_t numchans=12;  
    const Int_t binmax=2000;
    if(chan>=numchans){
      cout<<"\nERROR: requested channel is higher than total channels available! \n\n";
      break;
    }
    TH1F *hadc=new TH1F(name,source,binmax,0,binmax);
    FILE *infile = fopen(source,"r");
    cout<<"\n\tREADING IN FILE ---> "<<source<<"\n\n";
    Int_t adc[numchans];
    char cevent[256], cmodule[256];
    Int_t evnum, ncols;
    while(1){
      ncols = fscanf(infile,"\n%5c %i", &cevent,&evnum);
      if(ncols < 0) break;
      fscanf(infile,"\n%4c", &cmodule);
      for(Int_t j=0; j<12; j++){
	fscanf(infile," %i", &adc[j]);
      }
      if(adc[chan]<=binmax) hadc->Fill(adc[chan]-ped);
      //hadc->Fill(adc[chan]-ped);
    }
    cout<<"\n\tDONE READING FILE ---> "<<source<<"\n\n";
    fclose(infile);  
  }
  
  // USER DEFINED
  if(daq==2){
    const Int_t numchans=tnumchans;  
    const Int_t binmax=tbinmax;
    if(chan>=numchans){
      cout<<"\nERROR: requested channel is higher than total channels available! \n\n";
      break;
    }
    TH1F *hadc=new TH1F(name,source,binmax,0,binmax);
    ifstream infile;
    infile.open(source);
    if(infile.is_open()){
      infile.clear();
      infile.seekg(0,ios::beg);
      cout<<"\n\tREADING IN FILE ---> "<<source<<"\n";
      Int_t adc[numchans];
      Int_t j=0, events;  
      Float_t n=0;	      
      while(infile>>n){
	if(j==numchans) j=0; 
	if(j==chan) hadc->Fill(n-ped); 
	j++;		
      }
      cout<<"\n\tDONE READING FILE ---> "<<source<<"\n\n";
      infile.close();  
    }
    else cout<<"\n********** Could Not Read in File --> "<<source<<endl;
  }
  
  
  
  return hadc;
}

//--- RETURN THE NUMBER OF COLUMNS IN A FILE
int getCol(char *datafile)
{
  int info=0;
  FILE *infile=fopen(datafile,"r");
  int temp=0,prev=0,first=1;
  while((temp=getc(infile))!=EOF){
    
    //std::cout<<dec<<temp<<"  ";
    
    // if the first character(s) is a space or tab ignore it
    if((temp==32 || temp==9) && first==1 && prev==0){
      first=0;
      prev=temp;
      continue;
    }
    first=0;
    // if there's a space or tab and -NOT- and previous space or tab
    if((temp==32 || temp==9) && prev!=32 && prev!=9){
      info++;
      prev=temp;
      continue;
    }
    // if there's a newline and -NOT- a previous space or tab
    if(temp==10 && prev!=32 && prev!=9){
      info++;
      break;
    }
    // if there's a newline and -YES- a previous space or tab
    if(temp==10 && (prev==32 || prev==9)) break;
    prev=temp;
  
  } // end of while loop
  
  fclose(infile);
  return info;
}

//--- RETURN THE NUMBER OF LINES IN A FILE
int getLines(char *datafile)
{
  int info=0;
  ifstream infile;
  infile.open(datafile);
  string line;
  while(getline(infile,line)) info++;
  infile.close();
  return info;
}


