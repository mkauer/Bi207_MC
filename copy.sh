#!/bin/bash

mydate='09.09.23'
dir="${mydate}_BiESAT"

if [ "$1" ];then
    mkdir -p ./results/$dir/$1
    cp -fv AA_*.* ./results/$dir/$1
else
    echo "specify a file name"
fi

