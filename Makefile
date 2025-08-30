CXX := g++
ACCFLAGS_1 := -O3  -march=native  -I eigen3/ -fopt-info-vec -fopenmp


run: main.cpp
	${CXX} ${ACCFLAGS_1} -o run main.cpp 


clean:
	rm -f *.o *.csv *.png run *.plt *.dat
