CXX = g++
CXXFLAGS = -Wall -Wextra
LIBS =	-lmpi -lpthread

SRC = $(wildcard *.cpp ../utils/*.cpp)
OBJ = $(SRC:.cpp=.o)
TARGET = main

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ)

deep_clean: clean
	rm -f $(TARGET)

