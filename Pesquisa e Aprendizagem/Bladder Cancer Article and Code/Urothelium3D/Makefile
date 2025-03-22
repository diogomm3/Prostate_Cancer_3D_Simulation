SHELL := /bin/bash
CXX = g++

DEPS = src/CPM3D.h src/disjoint_set.h src/mt19937ar.h
OBJ = obj/CPM3D.o obj/disjoint_set.o obj/cpm_uro3D.o obj/mt19937ar.o
obj/%.o: src/%.cpp $(DEPS)
	$(CXX) -c -g -O3 -march=native -o $@ $<

run: $(OBJ)
	$(CXX)  $(OBJ) -g -O3 -march=native -o cpm_uro3D.out

clean:
	rm -f cpm_uro3D.out $(OBJ)
