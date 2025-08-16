# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra

# External libraries
# LHAPDF_CXXFLAGS := $(shell lhapdf-config --cflags)
# LHAPDF_LIBS     := $(shell lhapdf-config --libs)
# ROOT_CXXFLAGS   := $(shell root-config --cflags)
# ROOT_LIBS       := $(shell root-config --libs)
GSL_LIBS        := -lgsl

# Sources and objects
SRC  := main.cpp vegas.cpp
OBJ  := $(SRC:.cpp=.o)

# Output executable
TARGET := program.exe

# Default rule
all: $(TARGET)

# Link
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(GSL_LIBS)

# Compile rules
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

vegas.o: vegas.cpp vegas.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJ) $(TARGET)

# Run the program
run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
