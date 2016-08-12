#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>

#include <sys/stat.h>
#include <sys/unistd.h>


bool isaDirectory( std::string );
bool isaFile ( std::string );
bool isaEPS ( std::string );
bool isMultiple( int, int );

/////////////////////////////////////////////////////
// find and add all EPS files

void findEPS ( std::ifstream * , std::ofstream * , const char * , const char * );

///////////////////////////////////////////////////
// text utilities

std::string removeUpperDir( const std::string & , const char *);

std::string replaceSlashes( const std::string & );

std::string replaceSpecialChar( const std::string & );

std::string extractFileName( const char * );

std::string removeExtension( const std::string & , const char *);
