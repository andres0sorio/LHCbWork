// $Id: PDFs.h,v 1.2 2006/11/21 16:16:19 cmclean Exp $
#ifndef PDFS_H 
#define PDFS_H 1

// Include files
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cmath>
#include <vector>

#include "Riostream.h"
#include <TROOT.h>
#include <TMath.h>

double TotalPDF(double *x, double *par);

#endif // PDFS_H
