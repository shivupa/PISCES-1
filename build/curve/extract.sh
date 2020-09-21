#! /bin/sh
echo "#dist ADC2 diag-ADC2"
for i in 2.5 3.5  4 4.2  4.5 5.5 6.5 7.5 8.5 9.5
do
modelpot=` grep "Hartree" "$i"ang/out.dat  | tail -1| awk '{print $(NF-1)}'`;
echo  "$i $modelpot "
done

