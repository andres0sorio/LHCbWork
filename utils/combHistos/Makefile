SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

################
##linux

CXX          = g++
RM           = rm -f 
AR	     = ld

##Flags

CXXFLAGS     = -O -Wall -fPIC

#CXXFLAGS     += -D_DEBUG

ARFLAGS       =

#############
# Directories
INCDIR       = ../inc
EXEDIR       = ../exe
LIBDIR       = ../lib

CLHEPDIR     =$(HOME)/hepsoft/clhep_1.9

##Dependencies

### Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += $(ROOTCFLAGS) -I$(CLHEPDIR)/include -I$(INCDIR)
CLHEPLIB     = -L$(CLHEPDIR)/lib -lCLHEP
endif

###MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
CXXFLAGS    += -I$(CLHEPDIR)/include -I$(ROOTSYS)/include -I$(INCDIR)
CLHEPLIB     = -L$(CLHEPDIR)/lib -lCLHEP-CC.1.8.2.0
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif

#-----

LIBS        += $(ROOTLIBS) $(CLHEPLIB)
GLIBS       += $(ROOTGLIBS)

TOOLSLIBS    = -L$(LIBDIR) -lTools

#############################

OBJS     = combHistograms.o

GARBAGE  = *.o

########################################################################

all:combHistograms

########################################################################

%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< 	

combHistograms : $(OBJS)
	$(CXX) -o $(EXEDIR)/$@ $(OBJS) $(LINKFLAGS) $(TOOLSLIBS) $(GLIBS) $(ROOTLIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(EXEDIR)/combHistograms




