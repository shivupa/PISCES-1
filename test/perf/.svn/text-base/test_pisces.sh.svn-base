#! /bin/bash

##this script would test if the current version produces correct results for various methods in pisces

export PISCES_HOME=/home/kjordan/vkv3/pisces/2013_oct
export  PISCES_EXE=$PISCES_HOME/build/pisces

#test files:  w45_pol3_sp.inp w4_pol3_opt_small_input.inp w4_sp.inp  w8_d2h.inp 
##
$PISCES_EXE $PISCES_HOME/test/w4_sp.inp > w4_sp.out
echo "File = w4_sp.inp; ebe check"
grep EBE w4_sp.out
echo "Reference:
       EBE = -203.783 meV "
##
$PISCES_EXE $PISCES_HOME/test/w8_d2h.inp > w8_d2h.out
echo "File = w8_d2h.inp; excited states check"
grep "State [0-9]" w8_d2h.out
echo " Reference:
State 0     -0.03584329 Hartree =     -975.35 meV
State 1     -0.01128936 Hartree =     -307.20 meV
State 2     -0.00021632 Hartree =       -5.89 meV "


$PISCES_EXE $PISCES_HOME/test/w6a_opt.inp > w6a_opt.out
echo "File = w6a_opt.inp; optimization check"
grep Etotal w6a_opt.out
echo "Reference:
Step1     Etotal   = -1984.4 meV
Step2     Etotal   = -1442.7 meV
Step3     Etotal   = -2003.75 meV
Step4     Etotal   = -2008.03 meV
Step5     Etotal   = -2010.8 meV
Step6     Etotal   = -2012.56 meV
Step7     Etotal   = -2013.37 meV
Step8     Etotal   = -2014.1 meV "






