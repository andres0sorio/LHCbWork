#include "procStats.h"

procStats::procStats() {
  stats = NULL;
}

procStats::procStats(AllStats *statsIn) {

  isSignal = false;
  isBackground = false;

  stats = statsIn;
  nSig = stats->signal.size();
  nBack = stats->backgs.size();
  
  if(nSig > 0) isSignal = true;
  else if(nBack > 0) isBackground = true;
  else {
    std::cout << "There must be an error on your file" << std::endl;
    exit(1);
  }
  
}

procStats::~procStats() {
}

void procStats::printStats() {
  std::cout << "procStats> Printing stats record" << std::endl;
  for(unsigned int i = 0; i < (stats->signal).size(); i++) std::cout << *(stats->signal[i]);
  for(unsigned int i = 0; i < (stats->backgs).size(); i++) std::cout << *(stats->backgs[i]);
}


void procStats::loadData() {
  
 
  if(isSignal) {
    for ( int k = 0; k < nSig; k++) {
      wwXsec.push_back(stats->signal[k]->nnWWXsec);
      zzXsec.push_back(stats->signal[k]->nnZZXsec);
      bkXsec.push_back(stats->signal[k]->backXsec);

      nnZZentries.push_back(stats->signal[k]->nnZZevt);
      nnWWentries.push_back(stats->signal[k]->nnWWevt);
      nnBKentries.push_back(stats->signal[k]->backevt);
      
      wwCPairingTrue.push_back(stats->signal[k]->trueCPWW);
      wwCPairingDet.push_back(stats->signal[k]->detCPWW);
      zzCPairingTrue.push_back(stats->signal[k]->trueCPZZ);
      zzCPairingDet.push_back(stats->signal[k]->detCPZZ);
      
      nEventsWWtrue.push_back(stats->signal[k]->trueIntMWW_WW);
      nEventsWWdet.push_back(stats->signal[k]->detIntMWW_WW);
      nEventsZZtrue.push_back(stats->signal[k]->trueIntMWW_ZZ);
      nEventsZZdet.push_back(stats->signal[k]->detIntMWW_ZZ);
      nEventsBKtrue.push_back(stats->signal[k]->trueIntMWW_BK);
      nEventsBKdet.push_back(stats->signal[k]->detIntMWW_BK);

    }
    
    std::cout << "procStats> Total signal files processed: " << nSig << std::endl;
  }
  
  if(isBackground) {
    for ( int k = 0; k < nBack; k++) {
      
      wwCPairingTrue.push_back(stats->backgs[k]->trueCPWW);
      wwCPairingDet.push_back(stats->backgs[k]->detCPWW);
      zzCPairingTrue.push_back(stats->backgs[k]->trueCPZZ);
      zzCPairingDet.push_back(stats->backgs[k]->detCPZZ);
      
      nEventsBKtrue.push_back(stats->backgs[k]->trueIntMWW_BK);
      nEventsBKdet.push_back(stats->backgs[k]->detIntMWW_BK);
    }
    
    std::cout << "procStats> Total backgrounds files processed: " << nBack << std::endl;
  }
  
}

void procStats::getResults() { 
  
  if(isSignal) {
    
    wwXsecAv = sumVector(wwXsec);
    zzXsecAv = sumVector(zzXsec);
    bkXsecAv = sumVector(bkXsec);

    zzEntriesTot = sumVector(nnZZentries);
    wwEntriesTot = sumVector(nnWWentries);
    bkEntriesTot = sumVector(nnBKentries);

    wwCPairingTrueAv = evalAvg(wwCPairingTrue);
    wwCPairingDetAv = evalAvg(wwCPairingDet);
    
    zzCPairingTrueAv = evalAvg(zzCPairingTrue);
    zzCPairingDetAv = evalAvg( zzCPairingDet);
    
    nEventsWWtrueAv = sumVector(nEventsWWtrue);
    nEventsWWdetAv = sumVector(nEventsWWdet);
    
    nEventsZZtrueAv = sumVector(nEventsZZtrue);
    nEventsZZdetAv = sumVector(nEventsZZdet);
    
    nEventsBKtrueAv = sumVector(nEventsBKtrue);
    nEventsBKdetAv = sumVector(nEventsBKdet);

  }

  else if(isBackground) {
    
    wwCPairingTrueAv = evalAvg(wwCPairingTrue);
    wwCPairingDetAv = evalAvg(wwCPairingDet);
    zzCPairingTrueAv = evalAvg(zzCPairingTrue);
    zzCPairingDetAv = evalAvg(zzCPairingDet);

    nEventsBKtrueAv = evalAvg(nEventsBKtrue);
    nEventsBKdetAv = evalAvg(nEventsBKdet);

  }


}

double procStats::sumVector(const std::vector<double> & avec) {

  double result;

  std::vector<double>::const_iterator itr;
  result = 0.0;

  for(itr = avec.begin(); itr != avec.end(); ++itr) {
    result += (*itr);
  }
  
  return result;


}

int procStats::sumVector(const std::vector<int> & avec) {

  int result;

  std::vector<int>::const_iterator itr;
  result = 0;

  for(itr = avec.begin(); itr != avec.end(); ++itr) {
    result += (*itr);
  }
  
  return result;

}


double procStats::evalAvg(const std::vector<double> & avec) {
  
  double result;
  
  std::vector<double>::const_iterator itr;
  result = 0.0;
  
  for(itr = avec.begin(); itr != avec.end(); ++itr) {
    result += (*itr);
  }
  
  result = result / (avec.size());

  return result;

}

void procStats::produceTable() {

  std::cout << "------------------------------------------------------------" << std::endl;
  
  if(isSignal) {
    
    std::cout << "\t \t WW \t \t ZZ \t \t BK" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "nnWW X sec \t " << wwXsecAv
	      << std::endl;
    std::cout << "nnZZ X sec \t " << "\t \t" << zzXsecAv
	      << std::endl;
    std::cout << "6fBK X sec \t " << "\t \t \t \t" << bkXsecAv
	      << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Correct pairing rates" << std::endl;
    std::cout << "Detector data \t " << wwCPairingDetAv
	      << std::endl;
    std::cout << "Detector data \t " << "\t \t" << zzCPairingDetAv
	      << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "M_WW integral - total events after selection" << std::endl;
    std::cout << "Detector data \t " << nEventsWWdetAv
	      << std::endl;
    std::cout << "Detector data \t " << "\t \t" << nEventsZZdetAv
	      << std::endl;
    std::cout << "Detector data \t " << "\t \t \t \t" << nEventsBKdetAv
	      << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Total entries after selection" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

  }

  else if(isBackground) {

    std::cout << "\t \t WW \t \t ZZ \t \t BK" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Correct pairing rates" << std::endl;
    std::cout << "Detector data \t " << wwCPairingDetAv
	      << std::endl;
    std::cout << "Detector data \t " << "\t \t" << zzCPairingDetAv
	      << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "M_WW integral - total events after selection" << std::endl;
    std::cout << "Detector data \t " << "\t \t \t \t" << nEventsBKdetAv
	      << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Total entries after selection" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

  }

}
