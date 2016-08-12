// $Id: PlotUtility.h,v 1.5 2006/11/24 22:12:55 aosorio Exp $
#ifndef PLOTUTILITY_H 
#define PLOTUTILITY_H 1

// Include files

#include "Utilities.h"
#include "RandIntegrate.h"
#include "PDFs.h"
#include "PDFsPlot.h"
#include "PDFsTests.h"
#include "EvenOdd_PDFs.h"
#include "ResModels.h"
#include "RandomDists.h"

/** @PlotUtility PlotUtility.h
 *  
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-07-25
 */

void PlotUtility(double , double );

double projectTo1D(double *, double * );

void PlotPDFs();

void checkResModel();

//void testingToy( int );

#endif // PRODUCEDATA_H
