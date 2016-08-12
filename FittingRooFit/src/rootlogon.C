void rootlogon() {

  //gSystem->AddIncludePath(" -I/afs/cern.ch/sw/lcg/external/GSL/1.8/$CMTCONFIG/include ");
  //gSystem->AddtLinkedLibs(" -L/afs/cern.ch/sw/lcg/external/GSL/1.8/$CMTCONFIG/lib -lgsl");
  gSystem->Load("libMathMore.so");
  

  //don't compile
  //   gSystem->Load("Utilities_C.so");
  //   gSystem->Load("Integrate_C.so");
  //   gSystem->Load("ResModels_C.so");
  //   gSystem->Load("RandIntegrate_C.so");
  //   gSystem->Load("PDFs_C.so");
  //   gSystem->Load("RandomDists_C.so");
  //   gSystem->Load("UMLlh_Fit_C.so");

  
  //// compile
  gROOT->LoadMacro("Utilities.C++");
  gROOT->LoadMacro("FitResults.C++");
  gROOT->LoadMacro("GetResults.C++");
  

}
