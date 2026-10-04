CXX      ?= g++
CXXFLAGS ?= -O3 -ffast-math -march=native -std=c++17 -Wall -Wextra -pthread
all: dsr_encoder dsr_decoder
dsr_encoder: dsr_encoder.cpp dsr.hpp dsr_dsp.hpp dsr_io.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
dsr_decoder: dsr_decoder.cpp dsr.hpp dsr_dsp.hpp dsr_io.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
clean:
	rm -f dsr_encoder dsr_decoder
