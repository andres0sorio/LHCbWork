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

GARBAGE      = *.o $(EXEDIR)/generate $(EXEDIR)/generateWC

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

OBJSTHREE   = Utilities.o PDFsBkg.o ResModels.o \
	RandomDists.o ConvolGenerator.o produceData.o generate.o

OBJSFOUR   = Utilities.o PDFsBkg.o ResModels.o RandomDists.o \
	ConvolGenerator.o produceData.o generateWC.o

all : generate generateWC

%.o : %.C 
	$(CXX) $(CXXFLAGS) -c $< 	

generate: $(OBJSTHREE)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSTHREE) $(LINKFLAGS) $(LIBS)

generateWC: $(OBJSFOUR)
	 $(CXX) -o $(EXEDIR)/$@ $(OBJSFOUR) $(LINKFLAGS) $(LIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(SOGARBAGE)
