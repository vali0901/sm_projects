CXX = mpicxx
CXXFLAGS = -Wall -Wextra
LIBS = -lmpi
INCLUDE = -I../utils

SRC = $(wildcard *.cpp ../utils/*.cpp)
OBJ = $(SRC:.cpp=.o)
TARGET = main

NNODES = 4

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(INCLUDE) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	mpiexec -n $(NNODES) $(TARGET)

clean:
	rm -f $(OBJ)

deep_clean: clean
	rm -f $(TARGET)

