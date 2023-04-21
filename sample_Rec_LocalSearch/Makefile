CCC = g++ -std=c++11 -fopenmp

all: clean build

clean:
	rm -f evaluateTwitter.o
	rm -f Rec_sample.o
	rm -f SampleData.o
	rm -f BFS.o
	rm -f fixedDouble.o
	rm -f kMSolution.o
	rm -f loadsave.o
	rm -f vectorRec_LocalSearch.o
	rm -f sample_RecLocalSearchEval

build: evaluateTwitter.o
	$(CCC) -fopenmp -o sample_RecLocalSearchEval main.cpp evaluateTwitter.o Rec_sample.o vectorRec_LocalSearch.o fixedDouble.o kMSolution.o loadsave.o SampleData.o BFS.o

evaluateTwitter.o: Rec_sample.o loadsave.o SampleData.o
	$(CCC) -fopenmp -c evaluateTwitter.cpp -o evaluateTwitter.o

vectorRec_LocalSearch.o: fixedDouble.o kMSolution.o
	$(CCC) -fopenmp -c vectorRec_LocalSearch.cpp -o vectorRec_LocalSearch.o

fixedDouble.o:
	$(CCC) -fopenmp -c fixedDouble.cpp -o fixedDouble.o

kMSolution.o:
	$(CCC) -fopenmp -c kMSolution.cpp -o kMSolution.o

loadsave.o:
	$(CCC) -fopenmp -c loadsave.cpp -o loadsave.o

Rec_sample.o: vectorRec_LocalSearch.o
	$(CCC) -c Rec_sample.cpp -o Rec_sample.o

SampleData.o: BFS.o
	$(CCC) -c SampleData.cpp -o SampleData.o

BFS.o:
	$(CCC) -c BFS.cpp -o BFS.o
