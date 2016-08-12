SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

GSLPATH      =/afs/cern.ch/sw/lcg/external/GSL/1.8/$(CMTCONFIG)
GSLINC       =$(GSLPATH)/include
GSLLIB       =$(GSLPATH)/lib

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

GARBAGE      = *.o $(EXEDIR)/fitData

#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS) -I$(GSLINC)
LINKFLAGS    =
endif

##Libraries
LIBS       += $(ROOTLIBS) -lMathMore -lMinuit

OBJSONE    = Utilities.o PDFs.o OneAngleFit.o fitData.o

all : fitData

%.o : %.C 
	$(CXX) $(CXXFLAGS) -c $< 	

fitData : $(OBJSONE)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSONE) $(LINKFLAGS) $(LIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)
