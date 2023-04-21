#!/bin/sh

reset=0
rflag=10
sampleBegin=1
sampleEnd=1
sampleStep=1
while getopts k:r:p:b:e:s: flag
do
    case "${flag}" in
        k) kflag=${OPTARG};;
        r) rflag=${OPTARG};;
        p) reset=1;;
        b) sampleBegin=${OPTARG};;
        e) sampleEnd=${OPTARG};;
        s) sampleStep=${OPTARG};;
        *);;
    esac
done

echo $reset
if [ $reset -eq 1 ]
then
	rm ./EvaluationResults/RecLocalSearchRandomSampling/*.txt
	rm ./EvaluationResults/kUFLLocalSearchRandomSampling/*.txt
	# A_r, C_r, B_r, A_num, C_num, B_num
	python3 main.py 1 2 100 5 100 95 2
	xdg-open sphere_plot.png
fi

until [ $sampleBegin -gt $sampleEnd ]
do
  startrec=1
  endrec=${rflag}
  until [ $startrec -gt $endrec ]
  do
    ../Rec_LocalSearch_random_sampling/RecLocalSearchEval_random_sampling ./ ${sampleBegin}
    startrec=$(( $startrec + 1 ))
  done

  if [ "${kflag}" != "" ] && [ ${kflag} -gt 0 ]
  then
    rm ./EvaluationResults/kUFLLocalSearchRandomSampling/k=5_lam=1.txt
    ../kUFL_LocalSearch_random_sampling/kUFL_LocalSearch_random_sampling ./ "${kflag}" "${sampleBegin}"
    mv ./EvaluationResults/kUFLLocalSearchRandomSampling/k=5_lam=1.txt ./EvaluationResults/kUFLLocalSearchRandomSampling/k=5_lam="$sampleBegin".txt
  fi
  sampleBegin=$(( $sampleBegin + $sampleStep ))
done



