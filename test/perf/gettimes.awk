#!/usr/bin/awk -f

function output(arr, filename) {
   printf "%d ", size >> filename ;
   for (i in arr) 
      printf "%e ", arr[i] >> filename ;
   printf "\n" >> filename ;
}


/^OpenMP: num_threads = / { size = $4; } 
/^This is the Input read from/ {
   if ( length(data) == 0 ) {
      sub(/\.inp:/, "", $NF);
      name = $NF;
   } else name = data;
}   

/time: ComputePotential/ { tpot[npot++] = $4; }
/time: DVR::Diagonalize/ { tdia[ndia++] = $4; }
/time: ComputeGradient/  { tgra[ngra++] = $4; }
/time: Total /           { total[0] = $4; }

END {
   output( tpot, name  ".pot.time" );
   output( tdia, name  ".gra.time" );
   output( tgra, name  ".dia.time" );
   output( total, name ".tot.time" );
}
