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
EXEDIR       = .
LIBDIR       = ../lib


##Dependencies

# Linux
ifeq "$(UNAME)" "Linux"
RANLIB       = ranlib
CXXFLAGS    += $(ROOTCFLAGS) -I$(INCDIR)
endif

#MAC OSX
ifeq "$(UNAME)" "Darwin"
RANLIB       = ranlib -s
#CXXFLAGS   += -Wunreachable-code
CXXFLAGS    += -I$(INCDIR) $(ROOTCFLAGS)
LINKFLAGS    = -Xlinker -flat_namespace -bind_at_load
endif


LIBS        += $(ROOTLIBS) $(CLHEPLIB)
GLIBS       += $(ROOTGLIBS)

OBJS     = putLabels.o main.o

GARBAGE  = *.o

########################################################################

all: putLabels

########################################################################

%.o : %.cc 
	$(CXX) $(CXXFLAGS) -c $< 	

putLabels : $(OBJS)
	$(CXX) -o $(EXEDIR)/$@ $(OBJS) $(LINKFLAGS) $(GLIBS) $(ROOTLIBS)

clean	:
	$(RM) $(GARBAGE)

cleanall:
	$(RM) $(GARBAGE) $(EXEDIR)/putLabels




