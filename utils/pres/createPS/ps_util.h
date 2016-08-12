#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>

#include <sys/stat.h>
#include <sys/unistd.h>


bool isaDirectory( std::string );
bool isaFile ( std::string );
bool isaEPS ( std::string );
bool isMultiple( int, int );

int findNextMultiple( int , int );
/////////////////////////////////////////////////////
// find and add all EPS files

int addEPS ( std::ifstream * , std::ofstream * , const char * );

void updateTeXFile ( std::ofstream *, int );
///////////////////////////////////////////////////
// text utilities

std::string removeUpperDir( const std::string & , const char *);

std::string replaceSlashes( const std::string & );

std::string replaceSpecialChar( const std::string & );

std::string extractFileName( const char * );


