CC			:= gcc
CXX			:= g++

MPICC			:= mpicc
MPICXX			:= mpicxx

CCFLAGS			:= -O3 -march=native -Wall -std=gnu11
CXXFLAGS		:= -O3 -march=native -Wall -std=c++0x
LDFLAGS			:= -lm -lpthread -fopenmp -lboost_iostreams

all: seq seq_fibonacci seq_bellmen MPI_syn bellmen

seq: SSSP_seq.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_fibonacci: SSSP_seq_fibonacci.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_bellmen: SSSP_bellmen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
bellmen: Optimize_bellmen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MPI_syn: MPI_syn.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
clean:
	rm -f seq seq_fibonacci seq_bellmen MPI_syn bellmen
