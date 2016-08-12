#ifndef STATSSTRUCT_H
#define STATSSTRUCT_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>


struct SignalStats {
  
  // io functions
  friend std::istream& operator>>(std::istream &istr, SignalStats &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, SignalStats &rhs);
  
  // generated record  
  std::string process;
  std::string option;
  
  int nnZZevt, nnWWevt, backevt;
  double nnZZXsec, nnWWXsec, backXsec;
  double trueCPZZ, trueCPWW, detCPZZ, detCPWW;
  double trueIntMWW_WW, detIntMWW_WW;
  double trueIntMWW_ZZ, detIntMWW_ZZ;
  double trueIntMWW_BK, detIntMWW_BK;

};

struct BackgStats {
  
  // io functions
  friend std::istream& operator>>(std::istream &istr, BackgStats &rhs);
  friend std::ostream& operator<<(std::ostream &ostr, BackgStats &rhs);
  
  // generated record  
  std::string process;
  std::string option;
  
  double trueCPZZ, trueCPWW, detCPZZ, detCPWW;
  double trueIntMWW_BK, detIntMWW_BK;
  
};

class AllStats {
public:
  AllStats() {
  }
  AllStats(AllStats &e) {
    signal = e.signal;
    backgs = e.backgs;
  }
  
  ~AllStats() {
    
    std::vector<SignalStats*>::iterator itr1 = signal.begin();
    
    while(itr1 != signal.end()) {
      delete *itr1;
      itr1++;
    }
    
    std::vector<BackgStats*>::iterator itr2 = backgs.begin();

    while(itr2 != backgs.end()) {
      delete *itr2;
      itr2++;
    }

  }
  
  std::vector<SignalStats *> signal;
  std::vector<BackgStats *> backgs;
  
};

struct statsLoader {
  std::istream *is;
  statsLoader(const char *fileName);
  ~statsLoader();
  AllStats *next_file();
};


#endif
