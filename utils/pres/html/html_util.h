#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>


bool isaDirectory( std::string );
bool isaFile ( std::string );
bool isaPicture ( std::string );

/////////////////////////////////////////////////////
// find and add all GIF files

void addPictures ( std::ifstream * , std::ofstream * , const char * );

///////////////////////////////////////////////////
// text utilities

std::string removeUpperDir( const std::string & , const char *);

std::string replaceSlashes( const std::string & );

std::string extractFileName( const char * );

std::string removeExtension( const std::string & , const char *);
