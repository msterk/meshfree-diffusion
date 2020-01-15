#!/bin/bash

# scan over parameter number_of_nodes (n), on a single core of host computer


#declare -a numNodesArray=(1000 1258 1584 1995 2511 3162 3981 5011 6309 7943 10000 12589 15848 19952 25118 31622 39810 50118 63095 79432 100000 125892 158489 199526 251188 316227 398107 501187 630957 794328 1000000);
#declare -a numNodesArray=(1000000 192721 12100 784 49);

numNodesArray=(1024 1156 1369 1521 1764 2025 2401 2704 3136 3600 4096 4761 5476 6241 7056 8281 9409 10816 12321 14161 16384 18769 21609 24964 28561 32761 37636 43264 49729 57121 65536 75076 86436 99225 114244 131044 150544 173056 198916 228484 262144 301401 345744 396900 456976 524176 602176 692224 793881 912025 1048576)

numRepeats=10
numCores=1
alphaSupport=13
irregularity=0.3
numNodes="\${numNodes}"
outputName="result_scan_1Core_${numNodes}"
# add echo as the first command when debugging this script
command='echo ./mfre ${outputName} n=${numNodes} l=$numRepeats i=${irregularity} s=${alphaSupport}'


# overrides to debug the script and mfree application
#
#
#numNodesArray=(1024 1156 1369 1521 1764 2025 2401 2704 3136)
#outputName="-"
#
#
# ###################################################



echo '# scaan over parameter ${numNodes} using the following command line:'
echo "${command}"
echo ""

for numNodes in "${numNodesArray[@]}"
do
    echo "processing $numNodes nodes" >&2
    echo -ne "$numNodes\t";
	eval ${command}
done

echo ""
echo "Done" >&2

