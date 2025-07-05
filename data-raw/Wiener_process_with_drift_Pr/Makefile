# Makefile for Wiener_process_with_drift_Pr.cpp

CXX = clang++
CXXFLAGS = -std=c++11 -O3 -march=native -flto
TARGET = Wiener_process_with_drift_Pr
SRC = Wiener_process_with_drift_Pr.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)

