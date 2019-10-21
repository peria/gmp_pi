CC=clang
CXX=clang++
CFLAGS=-O2
CXXFLAGS=$(CFLAGS) -std=c++14
LDFLAGS=-lgmp -lm

all: pi_c pi_cc

pi_c: chudnovsky.c Makefile
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

pi_cc: chudnovsky.cc Makefile
	$(CXX) $(CXXFLAGS) -o $@ $< -lgmpxx $(LDFLAGS)

