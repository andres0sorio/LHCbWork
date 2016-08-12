void rootlogon() {

//don't compile
//  gSystem->Load("Utilities_C.so");
//  gSystem->Load("Integrate_C.so");
//  gSystem->Load("ResModels_C.so");
//  gSystem->Load("RandIntegrate_C.so");
//   gSystem->Load("PDFs_C.so");
//   gSystem->Load("RandomDists_C.so");
  
//// compile
  gROOT->LoadMacro("Utilities.C++");
  gROOT->LoadMacro("Integrate.C++");
  gROOT->LoadMacro("ResModels.C++");
  gROOT->LoadMacro("RandIntegrate.C++");
  gROOT->LoadMacro("Amplitudes.C++");
  gROOT->LoadMacro("PDFs.C++");
  gROOT->LoadMacro("RandomDists.C++");
  gROOT->LoadMacro("TestUtility.C++");

}
