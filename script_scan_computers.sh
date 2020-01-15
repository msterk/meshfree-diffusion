#!/bin/bash

# 36 computers, each has 4 cores; hosts contains this information in an open-mpi format
# scan over computers 1, 2, ... 36
# scan over parameter number_of_nodes (n): [100,400,1600,6400,10000, 22500,40000]


#declare -a numNodesArray=(1000 1258 1584 1995 2511 3162 3981 5011 6309 7943 10000 12589 15848 19952 25118 31622 39810 50118 63095 79432 100000 125892 158489 199526 251188 316227 398107 501187 630957 794328 1000000);

declare -a numNodesArray=(1000000 192721 12100 784 49);
numRepeats=10

#for (( numComputers=15; numComputers>0; numComputers-- ))
for (( numComputers=35; numComputers>1; numComputers-- ))
do
	numCores=$numComputers
    # radius is to be calculated from the number of nodes
    radius=2.5
	echo $numComputers
    for numNodes in "${numNodesArray[@]}"
    do
        echo -ne $numNodes"\t";
        mpirun -hostfile hosts_1core -n $numCores -bind-to-core ./mfre - n=$numNodes z=1 l=$numRepeats s=$radius
    done
done
echo 1
for numNodes in "${numNodesArray[@]}"
do
	echo -ne $numNodes"\t";
	./mfre - n=$numNodes z=0 l=$numRepeats s=$radius
done
