BOOST = # path to boost installation direction, need to contain both an "include" direction and a "lib" directory
EIGEN = # path to eigen library, need to contain a directory named "Eigen"
JELLYFISH = # path to jellyfish library, need to contain a directory named "jellyfish"
Spline = # path to the directory containing spline.h
Pbar = # path to progress-cpp directory, need to contain an "include" directory

CC = gcc
CXX = g++
CXXFLAGS = -std=c++11 $(INCLUDES) -g -fopenmp -O2
INCLUDES = -I $(BOOST) -I $(EIGEN) -I $(JELLYFISH)  -I $(Spline) -I $(Pbar)
LDADD = 
LDLIBS = -lz -lm -fopenmp -lboost_iostreams 

SRCS_init = utils/testgraph.cpp utils/Transcript.cpp utils/GeneGraph.cpp utils/GCbias.cpp utils/SEQbias.cpp utils/POSbias.cpp utils/CombinedBias.cpp
SRCS_bias = utils/pathbias.cpp utils/GCbias.cpp utils/SEQbias.cpp utils/POSbias.cpp utils/CombinedBias.cpp utils/Transcript.cpp utils/GeneGraph.cpp

all: bin/testgraph bin/pathbias

bin/testgraph: $(subst .cpp,.o,$(SRCS_init))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

bin/pathbias: $(subst .cpp,.o,$(SRCS_bias))
	mkdir -p bin
	$(CXX) -o $@ $^ $(LDADD) $(LDLIBS)

clean:
	rm -f bin/testgraph utils/*.o
