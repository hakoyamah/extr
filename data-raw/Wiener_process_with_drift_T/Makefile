CXX = clang++
CXXFLAGS = -std=c++11 -O3 -march=native -flto

SRC1 = Wiener_process_with_drift_T.cpp

TARGETS = Wiener_process_with_drift_T

all: $(TARGETS)

lande_orzack_cpp: $(SRC1)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC1)

clean:
	rm -f $(TARGETS)

