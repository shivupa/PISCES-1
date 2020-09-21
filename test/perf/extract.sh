#!/bin/bash

if [ $# == 0 ] ; then
    read -p "Use all .op files? (y/n) " x
    if [ "$x" != 'y' ] ; then exit 0; fi
    oplist=*.op
    read -p "Data name? (blank for default) " name
else 
    oplist=$@
fi

for x in $oplist ; do 
    echo $x
    gettimes.awk  -v data="$name" $x
done
