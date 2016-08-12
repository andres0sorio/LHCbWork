# Programs
CXX          = g++
CXXFLAGS     = -O1 -g -pg
#CXXFLAGS     = -O -Wall -fPIC
#CXXFLAGS     = -D_DEBUG


RM           = rm -f 

# Assign or Add variables
EXEDIR       = ../../exe/rootfiles
#EXEDIR        = .

PSOBJ     = ps_util.o createPS.o

GARBAGE      = *.o $(EXEDIR)/createPS

all : createPS

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< 	

createPS : $(PSOBJ)
	$(CXX) -o $(EXEDIR)/$@ $(PSOBJ)
	cp head_tex.dat $(EXEDIR)/.

clean	:
	$(RM) $(GARBAGE)





