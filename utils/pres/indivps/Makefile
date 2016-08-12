# Programs
CXX          = g++
CXXFLAGS     = -O1 -g -pg
#CXXFLAGS     = -O -Wall -fPIC
#CXXFLAGS     = -D_DEBUG


RM           = rm -f 

# Assign or Add variables
EXEDIR       = ../../exe/output
#EXEDIR        = ./

PSOBJ     = ps_util2.o createIndPS.o

GARBAGE      = *.o $(EXEDIR)/createIndPS

all : createIndPS

%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -c $< 	

createIndPS : $(PSOBJ)
	$(CXX) -o $(EXEDIR)/$@ $(PSOBJ)

clean	:
	$(RM) $(GARBAGE)





