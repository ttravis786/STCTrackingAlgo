CXX = g++
CXXFLAGS = -Iinclude -Wall -g `root-config --cflags`
LDFLAGS = `root-config --libs`
SRC = src/main.cpp
OBJ = src/main.o  # Specify the object file path
TARGET = my_program

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ `root-config --libs`

src/%.o: src/%.cpp  # Update the pattern to match src directory
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
