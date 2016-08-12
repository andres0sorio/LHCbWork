# Programs
CXX          = g++
CXXFLAGS     = -O1 -g -pg
#CXXFLAGS     = -O -Wall -fPIC
#CXXFLAGS     = -D_DEBUG


RM           = rm -f 

# Assign or Add variables
EXEDIR       = ../../exe/output
#EXEDIR        = ./

HTMLOBJ     = html_util.o createHTML.o

GARBAGE      = *.o $(EXEDIR)/createHTML

all : createHTML

%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< 	

createHTML : $(HTMLOBJ)
	$(CXX) -o $(EXEDIR)/$@ $(HTMLOBJ)

clean	:
	$(RM) $(GARBAGE)





