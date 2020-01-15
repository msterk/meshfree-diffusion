#!/bin/bash

# 36 computers, each has 4 cores; hosts contains this information in an open-mpi format
# scan over computers 1, 2, ... 36
# scan over parameter number_of_nodes (n): [100,400,1600,6400,10000, 22500,40000]


declare -a numNodesArray=(1000 1258 1584 1995 2511 3162 3981 5011 6309 7943 10000 12589 15848 19952 25118 31622 39810 50118 63095 79432 100000 125892 158489 199526 251188 316227 398107 501187 630957 794328 1000000);

	numComputers=1
    numCores=1
    # radius is to be calculated from the number of nodes
    radius=2.5
    echo "experimenting with "$numComputers" computers, "$numCores" cores"
    for numNodes in "${numNodesArray[@]}"
    do
        echo -ne $numNodes"\t";
        ./mfre - n=$numNodes l=10 s=$radius
    done
