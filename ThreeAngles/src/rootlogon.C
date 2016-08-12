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
  gROOT->LoadMacro("Integrate.C++");
  gROOT->LoadMacro("ResModels.C++");
  gROOT->LoadMacro("RandIntegrate.C++");
  gROOT->LoadMacro("PDFs.C++");
  gROOT->LoadMacro("EvenOdd_PDFs.C++");
  gROOT->LoadMacro("PDFsTests.C++");
  gROOT->LoadMacro("PDFsPlot.C++");
  gROOT->LoadMacro("UMLlh_Fit.C++");
  
  
  //gROOT->LoadMacro("RandomDists.C++"); // <- problem with headers
  //gROOT->LoadMacro("produceData.C++"); // <- depends on RandomDists.C
  
  //gROOT->LoadMacro("PlotUtility.C++"); // <- depends on RandomDists.C
  
  gROOT->LoadMacro("FitResults.C++");
  gROOT->LoadMacro("GetResults.C++");
  

}
