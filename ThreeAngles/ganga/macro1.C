void macro1(const char *in, const char *opt, const char *ccnevt) {

     gSystem->Load("Utilities_C.so");
     gSystem->Load("Integrate_C.so");
     gSystem->Load("ResModels_C.so");
     gSystem->Load("RandIntegrate_C.so");
     gSystem->Load("PDFs_C.so");
     gSystem->Load("EvenOdd_PDFs_C.so");
     gSystem->Load("PDFsTests_C.so");
     gSystem->Load("RandomDists_C.so");
     gSystem->Load("UMLlh_Fit_C.so");
     gSystem->Load("produceData_C.so");
     gSystem->Load("PlotUtility_C.so");

     int nevt        = atoi(ccnevt);

     produceData(in,opt,nevt);

}
