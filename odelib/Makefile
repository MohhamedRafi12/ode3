# build the ode solver library
# both shared and static libraries are built

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS   = $(ROOTCFLAGS) -I$(ODELIB)/inc -Wall -O3 -fPIC
LIBS       = $(ROOTLIBS)
GLIBS      = $(ROOTGLIBS)
GXX	   = g++ $(CXXFLAGS)



all: RKn.o libode.a libode.so

RKn.o: RKn.cpp RKn.hpp
	$(GXX) -c RKn.cpp

libode.a: RKn.o
	ar crs libode.a RKn.o $(ROOTLIBS)

libode.so: RKn.o
	g++ RKn.o -shared -o libode.so $(ROOTLIBS)


clean:
	rm -f RKn.o libode.a libode.so
