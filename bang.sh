#!/bin/bash

if [ "$1" ];then
    ./combineMC.sh $1
    ./run.sh $1
    ./copy.sh $1
    echo -e "\n\n DONE \n\n\n"
else
    echo "specify a file name"
fi

