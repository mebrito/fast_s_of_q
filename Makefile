# Compiler and flags
CXX=g++
CXXFLAGS=-std=c++23 -O3 -fopenmp
TARGET=fast_s_of_q
SOURCES=fast_s_of_q.cpp

# Default target
all: $(TARGET)

# Build the target
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

# Clean target
clean:
	rm -f $(TARGET)

# Rebuild target
rebuild: clean all

# Phony targets
.PHONY: all clean rebuild
