all: default

default: main.cpp bfs.cpp
	g++ -I../ -std=c++11 -fopenmp -O3 -g -o bfs main.cpp bfs.cpp ../common/graph.cpp ref_bfs.o
clean:
	rm -rf bfs  *~ *.*~
