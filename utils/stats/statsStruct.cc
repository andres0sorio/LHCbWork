#include "statsStruct.h"

////////////////////////////////////////////
//operator definitions
//////////////////////
//1
std::istream& operator>>(std::istream &istr, SignalStats &rhs) {
  istr >> rhs.process >> rhs.option
       >> rhs.nnZZevt >> rhs.nnZZXsec
       >> rhs.nnWWevt >> rhs.nnWWXsec
       >> rhs.backevt >> rhs.backXsec
       >> rhs.trueCPZZ  >> rhs.trueCPWW >> rhs.detCPZZ >> rhs.detCPWW
       >> rhs.trueIntMWW_WW >> rhs.detIntMWW_WW
       >> rhs.trueIntMWW_ZZ >> rhs.detIntMWW_ZZ
       >> rhs.trueIntMWW_BK >> rhs.detIntMWW_BK;

  return istr;
}

std::ostream& operator<<(std::ostream &ostr, SignalStats &rhs) {
  ostr << rhs.process << " " << rhs.option
       << " " << rhs.nnZZevt << " " << rhs.nnWWevt << " " << rhs.backevt
       << " " << rhs.nnZZXsec << " " << rhs.nnWWXsec << " " << rhs.backXsec
       << " " << rhs.trueCPZZ  << " " << rhs.trueCPWW << " " << rhs.detCPZZ << " " << rhs.detCPWW
       << " " << rhs.trueIntMWW_WW << " " << rhs.detIntMWW_WW
       << " " << rhs.trueIntMWW_ZZ << " " << rhs.detIntMWW_ZZ
       << " " << rhs.trueIntMWW_BK << " " << rhs.detIntMWW_BK
       << std::endl;
  return ostr;
}

//2

std::istream& operator>>(std::istream &istr, BackgStats &rhs) {
  istr >> rhs.process >> rhs.option
       >> rhs.trueCPZZ  >> rhs.trueCPWW >> rhs.detCPZZ >> rhs.detCPWW
       >> rhs.trueIntMWW_BK >> rhs.detIntMWW_BK;
  
 return istr;
}

std::ostream& operator<<(std::ostream &ostr, BackgStats &rhs) {
  ostr << rhs.process << " " << rhs.option
       << " " << rhs.trueCPZZ  << " " << rhs.trueCPWW << " " << rhs.detCPZZ << " " << rhs.detCPWW
       << " " << rhs.trueIntMWW_BK << " " << rhs.detIntMWW_BK
       << std::endl;
  return ostr;
}

///////////////////////////////////////////////
// statssLoader
///////////////////////////////////////////////
statsLoader::statsLoader(const char *fileName) {
  is = new std::ifstream(fileName);
  if(!is) {
    std::cout << "statsLoader> could not open input file> " << fileName << std::endl;
    exit(1);
  }
  
  //  char header[256]; // get rid off of the header - 
  //for(int i = 0 ; i<0; i++) { 
  //  is->getline(header,256);
  // }
  
}

statsLoader::~statsLoader() {
  delete is;
}

AllStats* statsLoader::next_file() {

  ////////////////////////////////////////////////////////////////////////////////////
  //  This was the most unkindest cut of all. -- William Shakespeare, "Julius Caesar"



  // create new statsEvent;
  AllStats *istat = new AllStats();
  
  ///////////////////////////////////////////////////////////
  // load in generated statistical results
  ////////////////////////////////////////
  int nFiles(0);
  int type(0);
  int i;
  char header[256];
  char c;

  if(!(*is)) 
    return NULL;

  for(int i = 0 ; i<1; i++) { 
    is->getline(header,256);
  }
  
  *is >> nFiles;
  *is >> type;
  
#ifdef _DEBUG
  std::cout << std::string(header) << std::endl;
  std::cout << "new file: " << nFiles << " " << type << std::endl;
#endif
  
  switch (type){
  case 1:
    {
      for(i=0; i<nFiles ; ++i) {
	// get rid off of the header - 
	is->getline(header,256);
	is->getline(header,256);
	is->getline(header,256);

	SignalStats *oneS = new SignalStats();
	*is >> *oneS;
	//std::cout << *oneS << std::endl;
	(istat->signal).push_back(oneS);
      }
      // std::cout << "statsLoader> Total files loaded : " << (istat->signal).size() << std::endl;
      *is >> c;
      break;
    }
  case 2:
    {
      for(i=0; i<nFiles ; i++) {
	// get rid off of the header - 
	is->getline(header,256);
	is->getline(header,256);
	is->getline(header,256);
	BackgStats *oneS = new BackgStats();
	*is >> *oneS;
	//std::cout << *oneS << std::endl;
	(istat->backgs).push_back(oneS);
      }
      //std::cout << "statsLoader> Total files loaded : " << (istat->backgs).size() << std::endl;
      *is >> c;
      break;
    }
  default:
    std::cout << "statsLoader> Type not  identified" << std::endl;
  }
  
  //////////////////////////////////////////////////////////
  
  return istat;
  
}




