# build the example code 
# to build your own program, add the EXE name in the EXES list below
# then include a target section following the example of one of the
# other programs below

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTFLAGS   = $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS) 
CXXFLAGS  += $(ROOTCFLAGS) -I$(ODELIB) -Wall -O3
LDFLAGS    = $(ROOTLIBS) $(ROOTGLIBS) -Wl,-rpath,$(ODELIB) -L$(ODELIB) -lode
GSLFLAGS   = -lgsl -lgslcblas
GXX	   = g++ $(CXXFLAGS)


EXES=baseball1 baseball2 RKnDemo RKnStep RKnTest projGSL vterm


all: $(EXES)

RKnDemo: RKnDemo.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

RKnStep: RKnStep.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

RKnTest: RKnTest.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

baseball1: baseball1.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS) 

baseball2: baseball2.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

vterm: vterm.cpp $(ODELIB)/RKn.hpp 
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS)

projGSL: projGSL.cpp
	$(GXX) $(CXXFLAGS) -o$@ $@.cpp $(LDFLAGS) $(GSLFLAGS)


clean:
	rm -f *.o *.so *.d *~ $(EXES)


