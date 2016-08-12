void rootlogon() 
{
  
  rootlogon(0);
  //don't compile: 0 compile: 1
  
}

void rootlogon(int option) {
	  
  gSystem->Load("libMathMore.so");
    
  //don't compile
  if (!option) {
    
  gSystem->Load("Utilities_C.so");
  gSystem->Load("ResModels_C.so");
  gSystem->Load("PDFsBkg_C.so");
  gSystem->Load("PDFs_C.so");
  gSystem->Load("PDFsWRes_C.so");
  gSystem->Load("Projections_C.so");

  } else {
    
    //// compile the basic stuff first
    gROOT->LoadMacro("Utilities.C++");
    gROOT->LoadMacro("ResModels.C++");
    gROOT->LoadMacro("PDFsBkg.C++");
    gROOT->LoadMacro("PDFs.C++");
    gROOT->LoadMacro("PDFsWRes.C++");

    //Only for plots
    gROOT->LoadMacro("Projections.C++");

  }
  
}
