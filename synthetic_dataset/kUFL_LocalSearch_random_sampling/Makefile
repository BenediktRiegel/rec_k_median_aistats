CC=g++ -std=c++14

all: clean build

clean:
	rm -f RandomSampling.o
	rm -f vectorEkUFL_LocalSearch.o
	rm -f evaluateTwitter.o
	rm -f kMSolution.o
	rm -f loadsave.o
	rm -f fixedDouble.o
	rm -f kUFL_LocalSearch.o
	rm -f vectorRec_LocalSearch.o
	rm -f kUFL_LocalSearch_random_sampling

build: evaluateTwitter.o kMSolution.o loadsave.o RandomSampling.o vectorEkUFL_LocalSearch.o fixedDouble.o vectorRec_LocalSearch.o
	$(CC) -fopenmp -o kUFL_LocalSearch_random_sampling main.cpp evaluateTwitter.o vectorEkUFL_LocalSearch.o RandomSampling.o kMSolution.o loadsave.o fixedDouble.o vectorRec_LocalSearch.o

RandomSampling.o:
	$(CC) -fopenmp -c RandomSampling.cpp -o RandomSampling.o

vectorEkUFL_LocalSearch.o:
	$(CC) -c vectorEkUFL_LocalSearch.cpp -o vectorEkUFL_LocalSearch.o

vectorRec_LocalSearch.o:
	$(CC) -c vectorRec_LocalSearch.cpp -o vectorRec_LocalSearch.o

evaluateTwitter.o:
	$(CC) -c evaluateTwitter.cpp -o evaluateTwitter.o

kMSolution.o:
	$(CC) -c kMSolution.cpp -o kMSolution.o

loadsave.o:
	$(CC) -c loadsave.cpp -o loadsave.o

fixedDouble.o:
	$(CC) -c fixedDouble.cpp -o fixedDouble.o
