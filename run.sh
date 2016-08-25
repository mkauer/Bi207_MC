#!/bin/bash


if [ "$1" ] && [ "$2" ];then
    echo "$1   $2"
    root -l -b -q "BiESAT.cxx(\"${1}\",${2})" 2>/dev/null | tee "AA_BiESAT-output.log"
elif [ "$1" ] && [ ! "$2" ];then
    echo "$1"
    root -l -b -q "BiESAT.cxx(\"${1}\",2)" 2>/dev/null | tee "AA_BiESAT-output.log"
else
    echo -e "\nERROR: specify a file name \n"
fi

