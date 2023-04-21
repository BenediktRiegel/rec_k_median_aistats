CC=g++ -std=c++14

all: clean build

clean:
	rm -f ReductionAlgo.o
	rm -f vectorEkUFL_LocalSearch.o
	rm -f evaluateTwitter.o
	rm -f kMSolution.o
	rm -f loadsave.o
	rm -f fixedDouble.o
	rm -f kUFL_LocalSearch

build: evaluateTwitter.o kMSolution.o loadsave.o ReductionAlgo.o vectorEkUFL_LocalSearch.o fixedDouble.o
	$(CC) -fopenmp -o kUFL_LocalSearch main.cpp evaluateTwitter.o vectorEkUFL_LocalSearch.o ReductionAlgo.o kMSolution.o loadsave.o fixedDouble.o

ReductionAlgo.o:
	$(CC) -fopenmp -c ReductionAlgo.cpp -o ReductionAlgo.o

vectorEkUFL_LocalSearch.o:
	$(CC) -c vectorEkUFL_LocalSearch.cpp -o vectorEkUFL_LocalSearch.o

evaluateTwitter.o:
	$(CC) -c evaluateTwitter.cpp -o evaluateTwitter.o

kMSolution.o:
	$(CC) -c kMSolution.cpp -o kMSolution.o

loadsave.o:
	$(CC) -c loadsave.cpp -o loadsave.o

fixedDouble.o:
	$(CC) -c fixedDouble.cpp -o fixedDouble.o
