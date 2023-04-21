SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

CCC = g++ -std=c++11 -fopenmp

CPLEX_STUDIO = /opt/ibm/ILOG/CPLEX_Studio201


CPLEXDIR      = $(CPLEX_STUDIO)/cplex
CONCERTDIR    = $(CPLEX_STUDIO)/concert

CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include


CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG
CCFLAGS = $(CCOPT) -I"$(CPLEXINCDIR)" -I"$(CONCERTINCDIR)"

CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
# Does windows have BINDIST?
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
# Warning BINDIST
CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) $(dynamic:yes=-L$(CPLEXBINDIR))

CPLEXLIB      = cplex$(dynamic:yes=2010)
CCLNFLAGS = -lconcert -lilocplex -l"$(CPLEXLIB)" -lm -lpthread -ldl




all: clean build

clean: 
	rm -f LP_Algo
	rm -f main.o
	rm -f evaluate.o
	rm -f kMSolution.o
	rm -f fixedDouble.o
	rm -f loadsave.o
	rm -f ReductionAlgo.o
	rm -f CharLiKUFL.o
	rm -f BFS.o
	rm -f bipartitroundingForCharikar.o


build:	main


# Warning BINDIST
main: main.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o LP_Algo main.o evaluate.o kMSolution.o loadsave.o ReductionAlgo.o CharLiKUFL.o BFS.o bipartitroundingForCharikar.o fixedDouble.o $(CCLNFLAGS)

main.o: main.cpp evaluate.o
	$(CCC) -c main.cpp -o main.o
	
evaluate.o: kMSolution.o loadsave.o fixedDouble.o ReductionAlgo.o
	$(CCC) -c evaluate.cpp -o evaluate.o

kMSolution.o:
	$(CCC) -c kMSolution.cpp -o kMSolution.o
	
loadsave.o:
	$(CCC) -c loadsave.cpp -o loadsave.o

fixedDouble.o:
	$(CCC) -c fixedDouble.cpp -o fixedDouble.o

ReductionAlgo.o: CharLiKUFL.o
	$(CCC) -c ReductionAlgo.cpp -o ReductionAlgo.o

CharLiKUFL.o: CharLiKUFL.cpp BFS.o bipartitroundingForCharikar.o fixedDouble.o
	$(CCC) -c $(CCFLAGS) CharLiKUFL.cpp -o CharLiKUFL.o

BFS.o:
	$(CCC) -c BFS.cpp -o BFS.o

bipartitroundingForCharikar.o:
	$(CCC) -c bipartitroundingForCharikar.cpp -o bipartitroundingForCharikar.o
