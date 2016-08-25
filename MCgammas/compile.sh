#!/bin/bash

function _res
{
    echo -e "\n\n==============================================="
    cat proba.car | grep -P res.*scriptme | awk \
	'{printf("\tRESOLUTION ==>  %s %\n",$3)}'
    echo -e "===============================================\n\n\n"
}

cernbin="${CERN}/${CERN_LEVEL}/bin"
nemolib="../deps/simul/prog/lib.Linux"

/bin/rm -fv proba.f proba.o proba.x

_res && sleep .8s

${cernbin}/nypatchy proba.car proba.f proba.cra .go

f77 -m32 -c proba.f

f77 -m32 -o proba.x proba.o -L${nemolib} -lgenbb \
`${cernbin}/cernlib_noshift geant,pawlib,graflib,\
packlib_noshift,mathlib,kernlib_noshift`


[ -f "proba.x" ] && ./proba.x && _res



