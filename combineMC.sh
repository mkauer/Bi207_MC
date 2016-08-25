#!/bin/bash

#root -l combineMC.cxx

if [ "$1" ];then
    root -l -b -q "_combineMC.cxx(\"${1}\")"  2>/dev/null
else
    root -l "_combineMC.cxx"  2>/dev/null
fi

