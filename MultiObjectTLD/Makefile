default: batchexample

all: batchexample sdlexample camexample

batchexample:
		g++ -Wall -Wno-write-strings -O3 -fopenmp batchExample.cpp -o batchExample

sdlexample:
		g++ -Wall -O3 -fopenmp `sdl-config --cflags` sdlExample.cpp -lSDL -lSDL_image -o sdlExample
		
camexample:
		g++ -Wall -O3 -fopenmp camExample.cpp `pkg-config opencv --cflags --libs` -o camExample

debug:
		g++ -Wall -Wno-write-strings -Wno-unknown-pragmas -g -pg batchExample.cpp -o batchExample

cleanall: clean cleanop

clean:
		rm -f sdlExample camExample batchExample

cleanop:
		rm -rf output/*

