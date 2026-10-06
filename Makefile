PROJECT := surface_field237
CXX := g++
CPPFLAGS += -Iinclude -Ithird_party/eigen3
CXXFLAGS += -std=c++20 -O2 -Wall -Wextra

SRC := src/topology.cpp src/validation.cpp src/qem.cpp src/simplify.cpp src/mesh_lod.cpp src/surface_field.cpp
OBJ := $(SRC:src/%.cpp=build/%.o)

.PHONY: all lib example test check clean

all: lib example test

lib: build/lib$(PROJECT).a

build:
	mkdir -p build bin

build/%.o: src/%.cpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

build/lib$(PROJECT).a: $(OBJ)
	$(AR) rcs $@ $(OBJ)

bin/open_surface: examples/open_surface.cpp build/lib$(PROJECT).a | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -Lbuild -l$(PROJECT) -o $@

bin/field_projection: examples/field_projection.cpp build/lib$(PROJECT).a | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -Lbuild -l$(PROJECT) -o $@

bin/test_main: tests/test_main.cpp build/lib$(PROJECT).a | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -Lbuild -l$(PROJECT) -o $@

example: bin/open_surface bin/field_projection
	./bin/open_surface
	./bin/field_projection

test: bin/test_main
	./bin/test_main

check: test example

clean:
	rm -rf build bin
