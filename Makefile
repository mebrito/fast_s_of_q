# Compiler and flags
CXX = g++
CXXSTD = -std=c++23
WARNFLAGS = -Wall -Wextra
OPTFLAGS = -O3
DEBUGFLAGS = -g
OMPFLAGS = -fopenmp
CXXFLAGS = $(CXXSTD) $(WARNFLAGS) $(OPTFLAGS) $(OMPFLAGS)

# Directories
SRCDIR = .
OBJDIR = obj
INCDIR = .

# Files
TARGET = fast_s_of_q
SOURCES = fast_s_of_q.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = s_of_q.hpp

# Default target
all: $(TARGET)

# Build the target
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
	
# Clean target
clean:
	rm -f $(TARGET) $(OBJECTS)

# Rebuild target
rebuild: clean all

# Phony targets
.PHONY: all clean rebuild
