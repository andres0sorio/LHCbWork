# Root variables
ROOTCFLAGS   =-L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     =-L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    =-L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

ECHO=echo
RM=rm -f	
CXX=g++

INCDIR=.
EXEDIR=../../exe
#EXEDIR=.

INCLUDES     =-I$(INCDIR) 
INCLUDES     +=-I$(ROOTSYS)/include

CXXFLAGS      =-O -Wall -fPIC
CXXFLAGS     +=$(ROOTCFLAGS)

#CXXFLAGS     +=-D_DEBUG
#CXXFLAGS    +=-D_NOASCII

LIBS         =$(ROOTLIBS)

OBJS         =main.o statsStruct.o procStats.o plotUtil.o

all : statsReader

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< 

statsReader : $(OBJS)
		$(CXX) $(OBJS) $(GLIBS) $(LIBS) -o $(EXEDIR)/$@

clean  : 
	$(RM) *.o
	$(RM) $(EXEDIR)/statsReader
