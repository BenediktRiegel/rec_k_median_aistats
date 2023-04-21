#!/bin/sh

reset=0
rflag=10
startd=2
endd=2
stepsize=1
while getopts k:r:s:e:t: flag
do
    case "${flag}" in
        k) kflag=${OPTARG};;
        r) rflag=${OPTARG};;
	s) startc=${OPTARG};;
	e) endc=${OPTARG};;
	t) stepsize=${OPTARG};;
        *);;
    esac
done

if [ $startc -lt 0 ]
then
	startc=0
fi

endc=$(( $endc + 1 ))
c=$startc
# rm -r ./num_clients_*
rm ./EvaluationResults/RecLocalSearch/*.txt
rm ./EvaluationResults/kUFLLocalSearch/*.txt
rm ./*.png
until [ $c -gt $endc ]
do
	# A_r, C_r, B_r, A_num, C_num, B_num, Dim
	python3 main.py 1 2 100 5 $c 20,20,20,20 2

	../Rec_LocalSearch/RecLocalSearchEval ./ ${rflag}
	mv ./EvaluationResults/RecLocalSearch/k=5_lam=1.txt ./EvaluationResults/RecLocalSearch/k=5_lam="$c".txt
	
	if [ $kflag -gt 0 ]
	then
		../kUFL_LocalSearch/kUFL_LocalSearch ./ "${kflag}"
		mv ./EvaluationResults/kUFLLocalSearch/k=5_lam=1.txt ./EvaluationResults/kUFLLocalSearch/k=5_lam="$c".txt
	fi
	echo "hi"
	c=$(( $c + $stepsize ))
done

# ./../plot.py meandiff minmaxdiff minmindiff mean meanratio

