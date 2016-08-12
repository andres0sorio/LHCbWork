#include "ps_util2.hpp"

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

bool haveSameName ( std::string objname , const char * aname ) {
  
  bool ans;
  
  if ( objname == std::string(aname)) ans = true;
  else ans = false;
  
  return ans;
  
}

std::string removeExtension( const std::string & filename, const char *ext) {

  std::string temp = filename;

  int iStart = filename.find(ext);  
  if(iStart < 0) iStart = 0;
  int iEnd = filename.length();
  return temp.erase(iStart,(iEnd-iStart));
  
}

bool isMultiple( int number, int factor) {
  
  bool ans;
  
  double r =  fmod(number+0.0,factor+0.0);

  if ( r == 0.0 ) ans = true;
  else ans = false;

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
    temp.replace(pos,1,"_");
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
    temp.replace(pos,1," ");
    iStart = pos;
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

void findEPS ( std::ifstream * currentdir , std::ofstream * targetfile , const char *option, const char * histoname ) {
  
  std::string object;
  bool dirtext = false;
  
  int counter = 1;
  
  while (currentdir->good())
    
    {
      
      *currentdir >> object;
      
      if(!currentdir->good()) break; //for some reason, the while does an extra loop at the end
      
      if(isaDirectory(object)) {
#if _DEBUG	
	std::cout << "dir found: " << object << std::endl;
#endif
	chdir(object.c_str());
	system("ls -p > list.txt");
  	std::ifstream * cdir;
	cdir = new std::ifstream("list.txt", std::ifstream::in);
	findEPS ( cdir , targetfile , option , histoname);
	delete cdir;
	system("rm list.txt");
	chdir("../");
	
      }
      
      else if(isaFile(object)) {
	
	if (isaEPS(object))  {

	  if(haveSameName(removeExtension(object,".eps") , histoname ) ) {
	    
	    const char *cwd;
	    unsigned int size = 256;
	    char buffer[size];
	    cwd = getcwd(buffer, size);
	    
	    std::string dirname = std::string(cwd);
	    
	    std::string sourcedir = removeUpperDir(dirname, option);
	    
	    std::string shortdirname = replaceSlashes(sourcedir);
	    
	    if(!dirtext) {
	      
	      //	      *targetfile << "\\clearpage"  << std::endl;
	      *targetfile << "\\section{" << replaceSpecialChar(shortdirname)  << "}" << std::endl;
	      dirtext = true;
	      
	    }
	    
	    if(isMultiple(counter,4)) *targetfile << "\\clearpage"  << std::endl;
	    
	    std::string source = dirname + std::string("/") + object;
	    
	    *targetfile << "\\begin{figure}[hp]" 
			<< std::endl
			<< "\\begin{center}" 
			<< std::endl
			<< "\\includegraphics[ scale=0.80]{"
			<< source
			<< "}"
			<< std::endl;
	    *targetfile << "\\caption{"
			<< replaceSpecialChar(object)
			<< "}" 
			<< std::endl;
	    *targetfile << "\\end{center}" << std::endl
			<< "\\end{figure}" << std::endl;
	    *targetfile << std::endl;
	    
	    //	  *targetfile << replaceSpecialChar(object) << std::endl;
	    
	    ++counter;
	    
	  }
	  
	}
	
      }
      
      else std::cout << "findEPS> What is this object?" << std::endl;
      
    }
  
  std::cout << "findEPS> Done." << std::endl;
  
}

