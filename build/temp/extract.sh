#! /bin/bash

echo "Pontential Diagonalize"
echo "C60"
for i in 1 2 4 8 16 
do 
   echo `grep  ComputePotential c60.out.$i | awk '{print $NF}'` `grep  Diagonalize c60.out.$i | awk '{print $NF}'`
done

echo "C240C60"
for i in 1 2 4 8 16
do 
   echo `grep  ComputePotential c240c60.out.$i | awk '{print $NF}'` `grep  Diagonalize c240c60.out.$i | awk '{print $NF}'`
done
