all: clean build

clean:
	rm -f evaluateTwitter.o
	rm -f fixedDouble.o
	rm -f kMSolution.o
	rm -f loadsave.o
	rm -f vectorRec_LocalSearch.o
	rm -f RecLocalSearchEval

build: evaluateTwitter.o
	g++ -fopenmp -o RecLocalSearchEval main.cpp evaluateTwitter.o vectorRec_LocalSearch.o fixedDouble.o kMSolution.o loadsave.o

evaluateTwitter.o: vectorRec_LocalSearch.o loadsave.o
	g++ -fopenmp -c evaluateTwitter.cpp -o evaluateTwitter.o

vectorRec_LocalSearch.o: fixedDouble.o kMSolution.o
	g++ -fopenmp -c vectorRec_LocalSearch.cpp -o vectorRec_LocalSearch.o

fixedDouble.o:
	g++ -fopenmp -c fixedDouble.cpp -o fixedDouble.o

kMSolution.o:
	g++ -fopenmp -c kMSolution.cpp -o kMSolution.o

loadsave.o:
	g++ -fopenmp -c loadsave.cpp -o loadsave.o