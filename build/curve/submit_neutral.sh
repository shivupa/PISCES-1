#! /bin/sh

for i in 2.5ang  3.5ang  4.2ang  4.5ang  4ang  5.5ang  6.5ang  7.5ang  8.5ang  9.5ang
do

cd $i/
echo "
BeginJob
 runtype = 1,
EndJob


BeginElectronPotential
  Polarization = 3,
  InternalParam = 1,
EndElectronPotential


BeginGridDef
 NoOfGridPoints = 80, 80, 80,
 Length = 81, 81, 81,
EndGridDef

BeginWaters" > input.dat
cat  ZMAT >> input.dat
echo "EndWaters" >> input.dat 
/ihome/kjordan/vkv3/pisces/latest/pisces/build/pisces input.dat > out.dat
cd ../

done

