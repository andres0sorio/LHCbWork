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
AR	     = ar cru

##Flags
CXXFLAGS     = -O0 -g -Wall -fPIC
#CXXFLAGS     += -D_DEBUG


INCDIR       = ../inc
SRCDIR       = ../src
LIBDIR       = ../lib


GARBAGE      = *.o

#################
##Dependencies
# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += -I$(CLHEP_BASE_DIR)/include -I$(INCDIR) $(ROOTCFLAGS)
CLHEPLIB     = -L$(CLHEP_BASE_DIR)/lib -lCLHEP
LINKFLAGS    =
endif

# MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
#CXXFLAGS   += -Wunreachable-code
CXXFLAGS    += -I$(CLHEPDIR)/include -I$(INCDIR) $(ROOTCFLAGS)
CLHEPLIB     = -L$(CLHEPDIR)/lib -lCLHEP
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif

##Libraries
LIBS       += $(ROOTLIBS) $(CLHEPLIB)


LIBOBJS    = string_utilities.o root_utilities.o \
	     addHistograms.o combHistograms.o printHistograms.o


all : libs


%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< 	

libs : $(LIBOBJS)
	$(AR) $(LIBDIR)/libTools.a $(LIBOBJS)
	$(RANLIB) $(LIBDIR)/libTools.a

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(LIBDIR)/libTools.a
