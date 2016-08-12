#ifndef PROCSTATS_H
#define PROCSTATS_H

#include <TROOT.h>
#include <TFile.h>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include "statsStruct.h"

class procStats {
  
 private:
  
  AllStats *stats;

  int nSig;
  int nBack;

 
  
 public:
  
  bool isSignal;
  bool isBackground;

  // define here all the variables you want to evaluate
  
  // Table Variables
  std::vector<double> wwXsec;
  std::vector<double> zzXsec;
  std::vector<double> bkXsec;

  std::vector<int> nnZZentries;
  std::vector<int> nnWWentries;
  std::vector<int> nnBKentries;

  std::vector<double> wwCPairingTrue;
  std::vector<double> wwCPairingDet;
  std::vector<double> zzCPairingTrue;
  std::vector<double> zzCPairingDet;
  
  std::vector<double> nEventsWWtrue;
  std::vector<double> nEventsWWdet;
  std::vector<double> nEventsZZtrue;
  std::vector<double> nEventsZZdet;
  std::vector<double> nEventsBKtrue;
  std::vector<double> nEventsBKdet;
  
  int zzEntriesTot;
  int wwEntriesTot;
  int bkEntriesTot;
  double wwXsecAv;
  double zzXsecAv;
  double bkXsecAv;
  double wwCPairingTrueAv;
  double wwCPairingDetAv;
  double zzCPairingTrueAv;
  double zzCPairingDetAv;
  double nEventsWWtrueAv;
  double nEventsWWdetAv;
  double nEventsZZtrueAv;
  double nEventsZZdetAv;
  double nEventsBKtrueAv;
  double nEventsBKdetAv;

   ///////////////////////////////////////////////////////
  // procCut constructors / destructor
  procStats();
  procStats(AllStats *);
  ~procStats();
  
  // Basic functions
  void printStats();
  
  void loadData();

  void getResults();

  double sumVector(const std::vector<double> &);
  
  int sumVector(const std::vector<int> &);

  double evalAvg(const std::vector<double> &);

  void produceTable();
  
};

#endif


