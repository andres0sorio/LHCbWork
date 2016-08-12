SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)
ifeq "$(CERN)" "/cern"
GSLPATH      =/afs/cern.ch/sw/lcg/external/GSL/1.8/$(CMTCONFIG)
GSLINC       =$(GSLPATH)/include
GSLLIB       =$(GSLPATH)/lib
else
GSLPATH      =$(LCG_release_area)/GSL/1.8/$(CMTCONFIG)
GSLINC       =$(GSLPATH)/include
GSLLIB       =$(GSLPATH)/lib
endif

################
##linux
CXX          = g++
RM           = rm -f 
AR	     = ar cru

##Flags
CXXFLAGS     = -O0 -g -Wall -fPIC
#CXXFLAGS     += -D_DEBUG


INCDIR       = .
SRCDIR       = .
LIBDIR       = .
EXEDIR       = ../exe

GARBAGE      = *.o $(EXEDIR)/stepOne $(EXEDIR)/stepTwo $(EXEDIR)/generate $(EXEDIR)/genparfiles

SOGARBAGE    = *.so

#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS) -I$(GSLINC)
LINKFLAGS    =
endif

# MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
#CXXFLAGS   += -Wunreachable-code
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS)
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif

##Libraries
LIBS       += $(ROOTLIBS) -lMinuit -lMathCore -lMathMore -L$(GSLLIB)

OBJSONE    = Utilities.o ResModels.o PDFs.o ProjectedPDFs.o DiffPDFs.o Distributions.o PDFsWacc.o \
        EvenOdd_PDFs.o Integrate.o RandomDists.o UMLlh_Fit.o stepOne.o

OBJSTWO    = Utilities.o ResModels.o PDFs.o ProjectedPDFs.o DiffPDFs.o Integrate.o \
	RandomDists.o UMLlh_Fit.o stepTwo.o

OBJSTHREE   = Utilities.o PDFs.o ProjectedPDFs.o DiffPDFs.o Integrate.o ResModels.o RandomDists.o \
	produceData.o generate.o

OBJSTEST    = Utilities.o genParFiles.o

OBJFITRES   = FitResults.o


all : stepOne stepTwo generate genparfiles fitresults

%.o : %.C 
	$(CXX) $(CXXFLAGS) -c $< 	

stepOne: $(OBJSONE)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSONE) $(LINKFLAGS) $(LIBS)

stepTwo: $(OBJSTWO)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSTWO) $(LINKFLAGS) $(LIBS)

generate: $(OBJSTHREE)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSTHREE) $(LINKFLAGS) $(LIBS)

genparfiles: $(OBJSTEST)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSTEST) $(LINKFLAGS) $(LIBS)

fitresults: $(OBJFITRES)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJFITRES) $(LINKFLAGS) $(LIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(SOGARBAGE)
