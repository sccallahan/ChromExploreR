#!/usr/bin/env bash
# store/run in 12chromHMM directory (or have folder with script in PATH)
# requires: 
# 	(1) perl script (2) hg19 genome (mouse or human) (3) cellmark file in directory with chromHMM output
num_states=15
num_loops=$((num_states+1))
# echo $num_loops
COUNTER=1
while [  $COUNTER -lt $num_loops ]; do
	echo The counter is $COUNTER
	perl chromStatesMatrixforClustering_particularStates.pl $num_states $COUNTER 10000 \
		namesChromXR.txt ./MYOUTPUT
	let COUNTER=COUNTER+1 
	done


# perl chromStatesMatrixforClustering_particularStates.pl 15 8 10000 sampleName.txt \
# /destination/directory