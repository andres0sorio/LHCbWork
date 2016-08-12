
# Programs
CXX          = g++
CXXFLAGS     = -O0 -Wall -g -fPIC -DGZSTREAM_NAMESPACE=gz

RM           = rm -f 

# Assign or Add variables
EXEDIR       = ../../exe
#EXEDIR       = .
INCDIR       = .

#gzstream 
GZSTREAM_DIR = $(HOME)/opt/gzstream

CXXFLAGS     += -I$(INCDIR) -I$(GZSTREAM_DIR)

GZSLIB       = -L$(GZSTREAM_DIR) -lgzstream
ZLIB         = -L/usr/lib -lz

LIBS         = $(GZSLIB) $(ZLIB)

CXXOBJS      = SimdetStruct.o SimdetEvent.o SimdetLoader.o main.o

GARBAGE      = *.o $(EXEDIR)/simdetAnal	

all : simdetAnal

%.o : %.cc 
	$(CXX) -c $< $(CXXFLAGS) 	

simdetAnal: $(CXXOBJS)
	$(CXX) -o $(EXEDIR)/$@ $(MAINOBJ) $(CXXOBJS) $(LIBS)

clean	:
	$(RM) $(GARBAGE)





