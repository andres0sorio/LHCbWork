void rootlogon() 
{
  
  rootlogon(1);
  //don't compile: 0 compile: 1
  
}

void rootlogon(int option) {
	  
  gSystem->Load("libMathMore.so");
  //gSystem->AddIncludePath("-I$HEPSOFT/lcg/external/GSL/1.8/$CMTCONFIG/include ");
  
  //gSystem->AddtLinkedLibs(" -L/afs/cern.ch/sw/lcg/external/GSL/1.8/$CMTCONFIG/lib -lgsl");
  
  //don't compile
  if (!option) {
    
  gSystem->Load("Utilities_C.so");
  gSystem->Load("ResModels_C.so");
  gSystem->Load("RandIntegrate_C.so");
  gSystem->Load("PDFs_C.so");
  gSystem->Load("ProjectedPDFs_C.so");
  gSystem->Load("DiffPDFs_C.so");
  gSystem->Load("EvenOdd_PDFs_C.so");
  gSystem->Load("PDFsWRes_C.so");
  gSystem->Load("RandomDists_C.so");
  gSystem->Load("ConvolGenerator_C.so");
  gSystem->Load("UMLlh_Fit_C.so");
  gSystem->Load("produceData_C.so");
  
  //gSystem->Load("PDFsTest_C.so");
  //gSystem->Load("PDFsPlot_C.so");

  } else {
    
    //// compile the basic stuff first
    gROOT->LoadMacro("Utilities.C++");
    gROOT->LoadMacro("ResModels.C++");
    gROOT->LoadMacro("RandIntegrate.C++");
    gROOT->LoadMacro("PDFs.C++");
    gROOT->LoadMacro("ProjectedPDFs.C++");
    gROOT->LoadMacro("DiffPDFs.C++");
    gROOT->LoadMacro("EvenOdd_PDFs.C++");
    gROOT->LoadMacro("PDFsWRes.C++");
    gROOT->LoadMacro("RandomDists.C++");
    gROOT->LoadMacro("ConvolGenerator.C++");
    gROOT->LoadMacro("UMLlh_Fit.C++");
    gROOT->LoadMacro("produceData.C++");
    
    //Only for testing
    //gROOT->LoadMacro("PDFsTests.C++");
    //gROOT->LoadMacro("PDFsPlot.C++");
    //has a problem with definitions - depends on Distributions
    //gROOT->LoadMacro("PlotUtility.C++"); 
  }
  
}
