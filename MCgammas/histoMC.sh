#!/bin/bash

if [ "$1" == "-q" ];then
    root -l -b -q "_histoMC.cxx"
elif [ "$1" ];then
    root -l "_histoMC.cxx(\"${1}\")"
else
    root -l '_histoMC.cxx'
fi

