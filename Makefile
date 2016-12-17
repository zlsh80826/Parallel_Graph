CC			:= gcc
CXX			:= g++

MPICC			:= mpicc
MPICXX			:= mpicxx

CCFLAGS			:= -O3 -march=native -Wall -std=gnu11
CXXFLAGS		:= -O3 -march=native -Wall -std=c++11
LDFLAGS			:= -lm -lpthread -fopenmp -lboost_iostreams

all: seq seq_fibonacci seq_bellmen syn bellmen synpr bellmenpc asyn

seq: SSSP_seq.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_fibonacci: SSSP_seq_fibonacci.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_bellmen: SSSP_bellmen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
bellmen: Optimize_bellmen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
bellmenpc: Optimize_bellmen_pc.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
syn: MPI_syn.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
synpr: MPI_syn_pr.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
asyn: MPI_asyn.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
clean:
	rm -f seq seq_fibonacci seq_bellmen syn bellmen synpr bellmenpc asyn
