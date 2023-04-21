#!/bin/sh

reset=0
rflag=10
while getopts k:r:p: flag
do
    case "${flag}" in
        k) kflag=${OPTARG};;
        r) rflag=${OPTARG};;
        p) reset=1;;
        *);;
    esac
done

echo $reset
if [ $reset -eq 1 ]
then
	rm ./EvaluationResults/RecLocalSearch/*.txt
	rm ./EvaluationResults/kUFLLocalSearch/*.txt
	# A_r, C_r, B_r, A_num, C_num, B_num, d
	python3 main.py 1 2 100 5 100 95 2
	# xdg-open sphere_plot.png
fi

../Rec_LocalSearch/RecLocalSearchEval ./ ${rflag}

if [ "${kflag}" != "" ] && [ ${kflag} -gt 0 ]
then
	../kUFL_LocalSearch/kUFL_LocalSearch ./ "${kflag}"
	rm ./*.png
	# ./../plot.py meandiff minmaxdiff minmindiff mean meanratio separatemeanratio
fi



