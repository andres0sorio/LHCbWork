# Root variables
ROOTCFLAGS   = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS     = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS    = -L$(ROOTSYS)/lib $(shell $(ROOTSYS)/bin/root-config --glibs)

# Programs
CXX          = g++
CXXFLAGS     = -O -Wall -fPIC -I$(ROOTSYS)/include
#CXXFLAGS     = -D_DEBUG

RM           = rm -f 

# Assign or Add variables
EXEDIR       = ../../exe/rootfiles

CXXFLAGS    += $(ROOTCFLAGS)

LIBS        += $(ROOTLIBS)
GLIBS       += $(ROOTGLIBS)

OBJECTS      = plot_util.o utilities.o printHistos.o

GARBAGE      = *.o $(EXEDIR)/printHistos

all : printHistos

%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -c $< 	

printHistos : $(OBJECTS)
	$(CXX) -o $(EXEDIR)/$@ $(OBJECTS) $(GLIBS) $(LIBS)

clean	:
	$(RM) $(GARBAGE)
