#include "ps_util.h"

bool isaDirectory ( std::string objname ) {
  
  bool ans;
  
  int check = objname.rfind('/');  
  if(check > 0) ans = true;
  else ans = false;
 
  return ans;

}

bool isaFile ( std::string objname ) {
  
  bool ans;
  
  int check1 = objname.rfind('/');  
  int check2 = objname.rfind('@');  
  if(check1 <= 0 && check2 <= 0 ) ans = true;
  else ans = false;
  
  return ans;

}

bool isaEPS ( std::string objname ) {
  
  bool ans;
  
  if (isaFile(objname) && !isaDirectory(objname)) {
    
    unsigned int loc = objname.find( ".eps", 0 );
    if( loc != std::string::npos ) ans = true;
    else ans = false;
    
  }
  
  else ans = false;
  
  return ans;
  
}

bool isMultiple( int number, int factor) {
  
  bool ans;
  
  double r =  fmod(number+0.0,factor+0.0);

  if ( r == 0.0 ) ans = true;
  else ans = false;

  return ans;
  
}


int findNextMultiple( int number, int factor) {
  
  int ans;
  
  double r =  fmod(number+0.0,factor+0.0);
  
  ans = 4 - int (r);
  
  return ans;
  
}

std::string removeUpperDir( const std::string & dirname, const char *option) {
  
  int iStart = dirname.find(option);  
  if(iStart < 0) iStart = 0;
  int iEnd = dirname.length();
  return dirname.substr(iStart,(iEnd-iStart));
  
}

std::string replaceSlashes( const std::string & name ) {
  
  std::string temp = name;
  
  int iStart = 0;
  int pos = name.find("/",iStart);
  
  while ( pos > 0 ) {
    temp.replace(pos,1," ");
    iStart = pos;
    pos = temp.find("/",iStart); 
  }
  
  return temp;
  
}

std::string replaceSpecialChar( const std::string & name ) {
  
  std::string temp = name;
  
  int iStart = 0;
  int pos = name.find("_",iStart);
  
  while ( pos > 0 ) {
    temp.insert(pos,"\\");
    iStart = pos+3;
    pos = temp.find("_",iStart); 
  }
  
  return temp;
}

std::string extractFileName( const char * name ) {
  
  std::string temp = std::string(name);
  
  int pos = temp.rfind('/');
  
  if(pos > 0) {
    temp.erase(pos);
  }
  
  pos = temp.rfind('/');
  
  if(pos > 0) {
    temp.erase(0,pos+1);
  }
  
  return temp;
  
}

// //////////////////////////////////////////////////////////
// //// Add pictures

int addEPS ( std::ifstream * currentdir , std::ofstream * targetfile , const char *option) {
  
  std::string object;
  bool dirtext = false;
  
  int counter = 0;
  int nCol = 1;

  //  std::vector<std::string> srcNames;
  
  while (currentdir->good())
    
    {
      
      *currentdir >> object;
      
      if(!currentdir->good()) break; //for some reason, the while does an extra loop at the end
      
      if(isaDirectory(object)) {
#ifdef _DEBUG	
	std::cout << "dir found: " << object << std::endl;
#endif
	chdir(object.c_str());
	system("ls -p > list.txt");
  	std::ifstream * cdir;
	cdir = new std::ifstream("list.txt", std::ifstream::in);
	int tot = addEPS ( cdir , targetfile , option );
	updateTeXFile(  targetfile , tot );
	delete cdir;
	system("rm list.txt");
	chdir("../");
	
      }
      
      else if(isaFile(object)) {
	
	if (isaEPS(object))  {
	  
#ifdef _DEBUG
	  std::cout << "eps found: " << object <<  std::endl;
#endif
	  
	  const char *cwd;
	  unsigned int size = 100;
	  char buffer[size];
	  cwd = getcwd(buffer, size);
	  
	  std::string dirname = std::string(cwd);

	  std::string sourcedir = removeUpperDir(dirname, option);
	  
	  std::string shortdirname = replaceSlashes(sourcedir);
	  
#ifdef _DEBUG
	  std::cout << "counter: " << counter << " ncol: " 
		    << nCol << std::endl;
#endif	  	  
	  
	  if(counter == 0 || isMultiple(counter,4)) {
	    nCol = 1;
	    
#ifdef _DEBUG
	    std::cout << "new Slide! " << std::endl; 
#endif

	    *targetfile << "\\lyxnewslide{wwsAnalysis}"  << std::endl;
	    *targetfile << replaceSpecialChar(shortdirname) << std::endl;
	    *targetfile << "\\begin{itemize}" << std::endl 
			<< "\\item \\textcolor{red}{wwsAnalysis results}"  << std::endl
			<< "\\end{itemize}" << std::endl
			<< "\\begin{tabular}{cc}" 
			<< std::endl;
	  }
	  
	  std::string source = dirname + std::string("/") + object;
	
	  //srcNames.push_back(source);
  
	  *targetfile << "\\includegraphics[scale=0.50]{"
		      << source
		      << "}";
	  
	  if( nCol == 2 ) {
	    
	    *targetfile << "\\tabularnewline" << std::endl;
	    nCol = 0;}
	  else if( nCol == 1) *targetfile << "&" << std::endl;
	  
	  ++counter;
	  
	  ++nCol;

	  if((counter == 4 || isMultiple(counter,4))) {
	    *targetfile << "\\end{tabular}" << std::endl
			<< std::endl;}
	  
	}
	
      }
      
      else std::cout << "addEPS> What is this object?" << std::endl;
      
    }
  
  std::cout << "addEPS> Done." << std::endl;

  return counter;
  
}

void updateTeXFile( std::ofstream * targetfile , int nFigs) {
  
  if ( isMultiple(nFigs,4) ) return;
  
  int nX;
  int nCol;
  
  nX = findNextMultiple ( nFigs, 4 );
  
  if(nX == 1) nCol = 2;
  else if (nX == 2) nCol = 1;
  else if (nX == 3) nCol = 2;
  else {}
  
  for ( int k = 0 ; k < nX; k++) { 
    
    *targetfile << "X "; 
    
    if( nCol == 2 ) {
      *targetfile << "\\tabularnewline" << std::endl;
      nCol = 0;}
    else if( nCol == 1) *targetfile << "&" << std::endl;
    ++nCol;
  }
  
  *targetfile << "\\end{tabular}" << std::endl
	      << std::endl;
  
  
}


