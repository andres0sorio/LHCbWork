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
EXEDIR       = .

GARBAGE      = *.o $(EXEDIR)/pdfValues

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
LIBS       += $(ROOTLIBS) -lMathMore -L$(GSLLIB)

OBJSONE    = ResModels.o Integrate.o PDFs.o Distributions.o \
	PDFsWacc.o PDFsTest.o pdfValues.o

all : pdfValues

%.o : %.C 
	$(CXX) $(CXXFLAGS) -c $< 	

pdfValues: $(OBJSONE)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSONE) $(LINKFLAGS) $(LIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE)
