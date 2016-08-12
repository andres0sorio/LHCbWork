#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "SimdetStruct.h"
#include "LeptonicOutputOrganizer.h"
#include "CutsLoader.h"
#include "ParamLoader.h"
#include "wwsLeptonEventProc.h"

#include "string_utilities.h"
#include "root_utilities.h"

int main(int iargv, const char **argv) {

  int arguments(0);

  arguments = iargv;
  
  if (arguments < 3 || arguments > 3) {
    std::cout << "usage: doLeptonicAnalysis [file] [nevents]" << std::endl;
    exit(1);
  }
  
  const char *fileName = argv[1];
  const int nEvent = atoi(argv[2]);
  
  std::cout << "wwsAnalysis> andres@hep.man.ac.uk" << std::endl;
  std::cout << "wwsAnalysis> reading file: " << fileName << std::endl;
  std::cout << "wwsAnalysis> events: " << nEvent << std::endl;
  
  std::string generator = getGeneratorName(fileName);
  
  std::string procname = getProcessName(fileName);
  
  std::string proctype = getProcessType(procname.c_str());
  
  proctype = std::string("signal");
  
  //////////////////////////////////////////////////
  //open input file
  SimdetLoader *sdl = new SimdetLoader(fileName);
  
  //open and read cuts from file
  CutsLoader *listofcuts = new CutsLoader("LeptonicSetofCuts.dat");
  Cuts *cuts = new Cuts();
  cuts->load_all_sets(listofcuts); 

  //open and read parameters from file
  ParamLoader *parameters = new ParamLoader("Parameters.dat");
  parameters->readParameters();

  ////////////////////////////
  // open output file
  TString outputfile = TString("./rootfiles/") + fileNameStub(fileName) + "_anal.root";
  TFile::TFile *outputFile = new TFile::TFile(outputfile,"recreate");  
  
  // define output data structure
  // first desginate the name of the structure - then name of the generator
  
  LeptonicOutputOrganizer *oc0 = new LeptonicOutputOrganizer("noCuts",
							     generator.c_str(),
							     proctype.c_str(),
							     procname.c_str());
  
  //LeptonicOutputOrganizer *oc1 = new LeptonicOutputOrganizer("withCuts", generator.c_str(), proctype.c_str(), procname.c_str());
  
  oc0->setParameters (parameters);
  //oc1->setParameters (parameters);
  
  SimdetEvent *event = sdl->next_event();

  int stage1(0),stage2(0),stage3(0);
  int stage4(0),stage5(0),stage6(0),stage7(0);
  int i(0);
  
  //while(event) {
  while(i < nEvent) {
    
#ifdef _DEBUG
    std::cout << "evt: " << i << " *" << std::endl;
#endif
    
    wwsLeptonEventProc *pevent = new wwsLeptonEventProc(event);
    
    pevent->setInternalParameters(800.00, generator, proctype);
    pevent->setExternalParameters(cuts,parameters);
    
    pevent->prepareEvent();

    ///////////////////////////////////
    //True MC simulation data analysis

    if(pevent->nTruePart > 0 ) {
      
      ++stage1;
      pevent->studyTrueData();
      pevent->applyProximityMethodTrue(0);
      oc0->fillWithTrueData(pevent);
      
    }
    
    ///////////////////////////////////
    //Detector simulation data analysis
    if(pevent-> nEflowObjs > 0 ) {
      
      ++stage2;
      pevent->studyDetectorData();
      //  if(pevent->isGoodForDetectorAnalysis()) {
      pevent->applyProximityMethod(0);
      oc0->fillWithDetectorData(pevent);
      //}
    }

#ifdef _DEBUG 
    // only if there is an interesting event
    //     if(pevent->interesting) {
    //       std::cout << "There is some peculiarity in this event" << std::endl;
    //       pevent->printEvent();
    //       delete pevent;
    //       delete event;
    //       //break;
    //     }
#endif			

    
    //delete event;
    delete pevent;
    
    
    event = sdl->next_event();
    
    ++i;
    
  }
  
  oc0->scaleHistograms();
  //oc1->importXsectionsFrom(oc0); // <<- This is exclusive for Whizard Data
  //oc1->scaleHistograms();
  oc0->evalDiffCrossSections();
  //oc1->evalDiffCrossSections();
  oc0->printStats("noCuts");
  //oc1->printStats("withCuts");
  
  // clean up
  
  outputFile->Write();
  outputFile->Close();
  
#ifdef _WARN  
  std::cout << "wwsAnalysis> Total events: " << i << std::endl;
  std::cout << "wwsAnalysis> Summary at stages: " 
	    << stage1 << " " 
	    << stage2 << " "
    	    << stage3 << " "
	    << stage4 << " "
	    << stage5 << " "
	    << stage6 << " "
	    << stage7 << " :: "
	    << stage4+stage6 << " "
	    << stage5+stage7 << std::endl;
  std::cout << "wwsAnalysis> eXterminated. " << std::endl;
#endif
  
  //delete cuts;
  //delete parameters;
  
  return 0;
  
}
