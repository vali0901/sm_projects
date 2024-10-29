CXX = g++
CXXFLAGS = -Wall -Wextra
LIBS = -lpthread
INCLUDE = -I../utils

SRC = $(wildcard *.cpp ../utils/*.cpp)
OBJ = $(SRC:.cpp=.o)
TARGET = main

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(INCLUDE) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ)

deep_clean: clean
	rm -f $(TARGET)
