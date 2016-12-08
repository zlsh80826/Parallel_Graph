CC			:= gcc
CXX			:= g++

MPICC			:= mpicc
MPICXX			:= mpicxx

CCFLAGS			:= -O3 -march=native -Wall -std=gnu11
CXXFLAGS		:= -O3 -march=native -Wall -std=c++0x
LDFLAGS			:= -lm

all: seq seq_fibonacci seq_set boost_test seq_bellmen

seq: SSSP_seq.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_fibonacci: SSSP_seq_fibonacci.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_set: SSSP_seq_set.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
boost_test: boost_test.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
seq_bellmen: SSSP_bellmen.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MS_MPI_dynamic_more: MS_MPI_dynamic_more.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MS_OpenMP_static: MS_OpenMP_static.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MS_OpenMP_dynamic: MS_OpenMP_dynamic.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MS_Hybrid_static: MS_Hybrid_static.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MS_Hybrid_dynamic: MS_Hybrid_dynamic.cc
	$(MPICXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
MS_seq: MS_seq.c
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $?
clean:
	rm -f seq seq_fibonacci seq_set boost_test seq_bellmen
