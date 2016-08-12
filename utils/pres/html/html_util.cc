#include "html_util.h"

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

bool isaPicture ( std::string objname ) {
  
  bool ans;
  
  if (isaFile(objname) && !isaDirectory(objname)) {
    
    unsigned int loc = objname.find( ".tiff", 0 );
    if( loc != std::string::npos ) ans = true;
    else ans = false;
    
  }
  
  else ans = false;
  
  return ans;
  
}

std::string removeUpperDir( const std::string & dirname, const char *option) {
  
  int iStart = dirname.find(option);  
  if(iStart < 0) iStart = 0;
  int iEnd = dirname.length();
  return dirname.substr(iStart,(iEnd-iStart));
  
}

std::string removeExtension( const std::string & filename, const char *ext) {

  std::string temp = filename;

  int iStart = filename.find(ext);  
  if(iStart < 0) iStart = 0;
  int iEnd = filename.length();
  return temp.erase(iStart,(iEnd-iStart));
  
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

void addPictures ( std::ifstream * currentdir , std::ofstream * targetfile , const char *option) {
  
  std::string object;
  bool dirtext = false;
  
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
	addPictures ( cdir , targetfile , option );
	delete cdir;
	system("rm list.txt");
	chdir("../");
	
      }
      
      else if(isaFile(object)) {
	
	if (isaPicture(object))  {
	  
#if _DEBUG
	  std::cout << "gif found: " << object <<  std::endl;
#endif
	  
	  const char *cwd;
	  unsigned int size = 100;
	  char buffer[size];
	  cwd = getcwd(buffer, size);
	  
	  std::string dirname = std::string(cwd);

	  std::string sourcedir = removeUpperDir(dirname, option);
	  
	  std::string shortdirname = replaceSlashes(sourcedir);
	  
	  if(!dirtext) {
	    
	    *targetfile << "<br>" << std::endl;
	    *targetfile << "<IMG src=\"directory.gif\" alt=\"diricon\" width=\"20\" height=\"17\">" << std::endl;
	    *targetfile << "<b>" << sourcedir << "</b><br>" << std::endl;
	    
	    dirtext = true;
	    
	  }
		 	  	  
	  std::string source = dirname + std::string("/") + object;

	  *targetfile << std::endl;
	  
	  *targetfile << "<A HREF=\""
		      << "javascript:openPopup ('"
		      << source << "','"
		      << object 
		      << "_"
		      << shortdirname
		      << "')\">"
		      << "&nbsp; &nbsp;"
		      << removeExtension(object,".tiff") 
		      << "</A>"
		      << "<br>" << std::endl;
	  
	}
	
      }
      
      else std::cout << "addPictures> What is this object?" << std::endl;
      
    }
  
  std::cout << "addPictures> Done." << std::endl;
  
}

