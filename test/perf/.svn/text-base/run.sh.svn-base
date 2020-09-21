#!/bin/bash

if [ $# == 0 ] ; then
    echo error: provide data file name
    exit 1
fi
data=$1
shift

if [ $# == 1 ] ; then
    echo info: using n `seq 1 8` 
    n=$(seq 1 8);
else 
    n=$*;
fi

rm -f $data.*.txt

for x in $n ; do 
    declare -x OMP_NUM_THREADS=$x
    ./pisces $data | tee $data.$x.op | awk '
function output(arr, filename) {
   printf "%d ", size >> filename ;
   for (i in arr) 
     printf "%e ", arr[i] >> filename ;
   printf "\n" >> filename ;
}
 /time: ComputePotential/ { tpot[npot++] = $4; }
 /time: DVR::Diagonalize/ { tdia[ndia++] = $4; }
 /time: ComputeGradient/  { tgra[ngra++] = $4; }
 /time: Total / { total[0] = $4; }
END {
   output( tpot, data  ".pot.time" );
   output( tdia, data  ".gra.time" );
   output( tgra, data  ".dia.time" );
   output( total, data ".tot.time" );
}' size=$x data=$data
done
