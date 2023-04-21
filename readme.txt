The synthetic dataset has its own readme.txt

We exclude the twitter dataset, due to its size.

First the algorithms have to be compiled. Use the their respective Makefiles, to do so. In kUFL_LocalSearch's case, the Makefile creates a new file called "kUFL_LocalSearch". 
This file is the algorithms executable.

To run an algorithm you need to provide it with two parameters:
1. The directory to the dataset
2. The number of runs (This is interesting for averaging the results)
An example could be:
./kUFL_LocalSearch/kUFL_LocalSearch ./data/congress 1

In this example, the results will be stored under:
./data/congress/EvaluationResults/kUFLLocalSearch

as .txt files, containing the used k and lambda in their names, e.g. "k=2_lam=0.3.txt"

The algorithms store their results in different directories. Here is a list:
kUFL_LocalSearch: 		EvaluationResults/kUFLLocalSearch/
kUFL_LP:			EvaluationResults/LP/
Rec_LocalSearch: 		EvaluationResults/RecLocalSearch/
sample_Rec_LocalSearch:		uniformSampleEvaluationResults/RecLocalSearch/
sparsified_kUFL_LocalSearch:	sparsifiedEvaluationResults/kUFLLocalSearch/

Make sure these directories exist within the data's directory, before starting the algorithm.


Datasets:

The data's directory should contain several files:
1. C.txt: 	contains one line with all the clients' ids, seperated by commas, e.g. 0,1,2,3,4
2. F.txt: 	contains one line with all the facilities' ids, seperated by commas, e.g. 1,4,5
3. dAtoC.txt:	the i'th line contains the distances to the first, second, third and so on client from the facility or client with id=i
4. k.txt:	contains one line with all the values for k, with which the algorithm should be executed, e.g. 2,3
5. lam.txt:	contains one line with all the values for lamda, with which the algorithm should be executed, e.g. 0,0.1,0.2


sample_Rec_LocalSearch and sparsified_kUFL_LocalSearch need the following additional files or file changes:
1. sample_k.txt:	analogous to k.txt
2. sample_lam.txt:	analogous to lam.txt
3. nearest_k.txt:	line i represents the i'th facility. Line i contains all facilities sorted by their distance to the i'th facility in an increasing order. Not containing the i'th facility.
4. nearest_f.txt:	contains one line with values separated by commas. The j'th value is the facility, closest to the j'th client.
5. sample_amounts:	contains one line with all the sample size, e.g. 6,8. sample_Rec_LocalSearch samples uniformly at random from C.