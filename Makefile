app: main.o
	g++ -o app main.o -lrpoly -llapack -lopenblas -lgfortran -pthread
main.o: main.cpp
	g++ -c main.cpp -o main.o -std=c++11

clean:
	rm app *.o
