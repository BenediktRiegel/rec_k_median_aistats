CCC = g++ -std=c++11 -fopenmp

all: clean build

clean:
	rm -f evaluateTwitter.o
	rm -f fixedDouble.o
	rm -f kMSolution.o
	rm -f loadsave.o
	rm -f vectorRec_LocalSearch.o
	rm -f RecLocalSearchEval_random_sampling

build: evaluateTwitter.o
	$(CCC) -o RecLocalSearchEval_random_sampling main.cpp evaluateTwitter.o vectorRec_LocalSearch.o fixedDouble.o kMSolution.o loadsave.o

evaluateTwitter.o: vectorRec_LocalSearch.o loadsave.o
	$(CCC) -c evaluateTwitter.cpp -o evaluateTwitter.o

vectorRec_LocalSearch.o: fixedDouble.o kMSolution.o
	$(CCC) -c vectorRec_LocalSearch.cpp -o vectorRec_LocalSearch.o

fixedDouble.o:
	$(CCC) -c fixedDouble.cpp -o fixedDouble.o

kMSolution.o:
	$(CCC) -c kMSolution.cpp -o kMSolution.o

loadsave.o:
	$(CCC) -c loadsave.cpp -o loadsave.o
