CXX = mpicxx
CXXFLAGS = -Wall -Wextra
LIBS = -lmpi

SRC = $(wildcard *.cpp ../utils/*.cpp)
OBJ = $(SRC:.cpp=.o)
TARGET = main

NNODES = 4

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	mpiexec -n $(NNODES) $(TARGET)

clean:
	rm -f $(OBJ)

deep_clean: clean
	rm -f $(TARGET)

