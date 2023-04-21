The directories "kUFL_LocalSearch", "kUFL_LocalSearch_random_sampling", "Rec_LocalSearch", "Rec_LocalSearch_random_sampling"
contain the source code for the algorithms.
Compile these first with the respective Makefile.

The directory "synthetic_cluster_sphere" contains:
k.txt: Define the k's used to compute a solution
lam.txt: Define the lambdas used to compute a solution
main.py: Creates the synthetic dataset
run.sh: Runs the algorithms located in "kUFL_LocalSearch" and "Rec_LocalSearch".
	The flag -k defines how often kUFL_LocalSearch should run
	The flag -r defines how often Rec_LocalSearch should run
	The flag -p (if set) deletes the EvaluationResults and computes a new dataset
./run.sh -r 10 -k 5 -p 1 -> runs Rec_LocalSearch 10 times, kUFL_LocalSearch 5 times and recomputes a new dataset

plot_clients.sh: Runs the algorithms located in "kUFL_LocalSearch" and "Rec_LocalSearch".
	The flags -k and -r are the same as above.
	The flag -s defines the min amount of clients that should be contained within the dataset.
	The flag -e defines the max amount of clients that should be contained within the dataset.
	The flag -t defines the stepsize for the client amount
./plot_clients.sh -r 10 -k 5 -s 10 -e 41 -t 10 -> runs Rec and kUFL local search with 10, 20, 30, 40 clients in the dataset.
Warning: Do not start with a flag -s 1
In the evaluation results "lam=" will know equal the amounts of clients in the dataset.

random_sampling.sh: Runs the random sampling versions of the two local search algorithms
	The flags -k, -r and -p are the same as above
	The flag b defines the min sample size.
	The flag e defines the max sample size.
	The flag s defines the step size
./random_sampling.sh -r 200 -k 200 -s 11 -e 32 -s 10 -> runs the algorithms each 200 times, with the sample size of 11, 21 and 31
Warning: Do not start with a flag -b 1
In the evaluation results "lam=" will know equal the sample size of the algorithms.


To change the dataset parameters you need to change the line containing "python3 main.py 1 2 100 5 100 95 2" to your desired values.
The order is as follows: python3 main.py good_r clients_r bad_r good_amount client_amount bad_amount dimensionality

plot.py is a simple python script to plot the results.
To plot the results run plot.py from within the directory of your dataset.
Examples: 
1. go into "synthetic_cluster_sphere" directory and run "./../plot.py mean" -> this creates a plot for the mean costs.
2. go into "synthetic_cluster_sphere" directory and run "./../plot.py random meanduration-s probability" 
-> this creates multiple plots for the mean runtime in seconds and the estimated probability of ending up with the optimum, based on the random sampling results.
For more options look into the method "def plot(ptypes)".