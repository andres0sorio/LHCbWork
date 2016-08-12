SHELL = /bin/sh
UNAME = $(shell uname)

# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs) 
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

CLHEP_BASE_DIR=/opt/cern/clhep/pro

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


##Dependencies

# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += $(ROOTCFLAGS) -I$(CLHEP_BASE_DIR)/include -I$(INCDIR)
CLHEPLIB     = -L$(CLHEP_BASE_DIR)/lib -lCLHEP
endif

#MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
#CXXFLAGS   += -Wunreachable-code
CXXFLAGS    += -I$(CLHEP_BASE_DIR)/include -I$(INCDIR) $(ROOTCFLAGS)
CLHEPLIB     = -L$(CLHEP_BASE_DIR)/lib -lCLHEP
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif


LIBS        += $(ROOTLIBS) $(CLHEPLIB)
GLIBS       += $(ROOTGLIBS)

TOOLSLIBS    = -L$(LIBDIR) -lTools


OBJS     = addHistograms.o

GARBAGE  = *.o

########################################################################

all:addHistograms

########################################################################

%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< 	

addHistograms : $(OBJS)
	$(CXX) -o $(EXEDIR)/$@ $(OBJS) $(LINKFLAGS) $(TOOLSLIBS) $(GLIBS) $(ROOTLIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(EXEDIR)/addHistograms

