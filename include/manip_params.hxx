/***********************************************************
FUNCTIONS FOR MANIPULATING THE PARAMETER NAMES
------------------------------------------------------------

VERSION: 08.11.02

change log: (most recent at top)
-----------------------------------
31oct + initialize function names

Matt Kauer (kauer@hep.ucl.ac.uk)
***********************************************************/


Int_t setParNames(TF1 *myfunc, string type){
  
  // set names for the bi207 fit parameters
  if(type=="bi207"){
    Int_t i=0;
    
    // expodential parValues
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"exp0");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"exp1");
    i++;
    
    // gamma 500
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Cn500");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Ce500");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Cs500");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Css500");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"CGn500");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"CGe500");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"CGs500");
    i++;
    
    // gamma 1000
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Cn1000");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Ce1000");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Cs1000");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Css1000");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"CGn1000");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"CGe1000");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"CGs1000");
    i++;
    
    // beta 482
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Gn482");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Ge482");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Gs482");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Gc482");
    i++;
    
    // beta 976
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Gn976");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Ge976");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Gs976");
    i++;
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"Gc976");
    i++;
    
    // offset
    if(i>MAXPARS) return i;
    myfunc->SetParName(i,"offset");
    
    return -1;
  }
}

