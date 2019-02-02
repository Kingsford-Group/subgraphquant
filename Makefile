BOOST = #/opt/local/stow/boost-1.55.0-gcc
EIGEN = #/home/congm1/ocean/oceancong02/Software/eigen-eigen-5a0156e40feb/
JELLYFISH = #/opt/local/stow/jellyfish-2.2.7/include/jellyfish-2.2.7/
Spline = #/home/congm1/ocean/oceancong02/Software/salmon-0.8.2/include/

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g
INCLUDES = -I $(BOOST) -I $(EIGEN) -I $(JELLYFISH)  -I $(Spline) 
LDADD = 
LDLIBS = -lz -lm -fopenmp -lboost_iostreams 

SRCS_init = utils/testgraph.cpp utils/Transcript.cpp utils/GeneGraph.cpp utils/GCbias.cpp utils/SEQbias.cpp utils/POSbias.cpp utils/CombinedBias.cpp
SRCS_bias = utils/biascorrection.cpp utils/GCbias.cpp utils/SEQbias.cpp utils/POSbias.cpp utils/CombinedBias.cpp

all: bin/testgraph bin/biascorrection

bin/testgraph: $(subst .cpp,.o,$(SRCS_init))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

bin/biascorrection: $(subst .cpp,.o,$(SRCS_bias))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

clean:
	rm -f bin/testgraph utils/*.o
