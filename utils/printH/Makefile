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

##Dependencies

# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += $(ROOTCFLAGS) -I$(CLHEPDIR)/include -I$(INCDIR)
CLHEPLIB     = -L$(CLHEPDIR)/lib -lCLHEP
endif

#MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
#CXXFLAGS   += -Wunreachable-code
CXXFLAGS    += -I$(CLHEPDIR)/include -I$(INCDIR) $(ROOTCFLAGS)
CLHEPLIB     = -L$(CLHEPDIR)/lib -lCLHEP
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif


LIBS        += $(ROOTLIBS) $(CLHEPLIB)
GLIBS       += $(ROOTGLIBS)

TOOLSLIBS    = -L$(LIBDIR) -lTools

OBJS     = printHistos.o

GARBAGE  = *.o $(EXEDIR)/printHistos

########################################################################

all:printHistos

########################################################################

%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< 	

printHistos : $(OBJS)
	$(CXX) -o $(EXEDIR)/$@ $(OBJS) $(LINKFLAGS) $(TOOLSLIBS) $(GLIBS) $(ROOTLIBS)
	cp $(EXEDIR)/printHistos $(HOME)/bin/.

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(EXEDIR)/printHistos




