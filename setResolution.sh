#!/bin/bash

function _err
{
    echo -e "\n ERROR: ${1} \n"
    exit 1
}

if [ "$1" ];then
    res="$1"
else
    res=6.5
fi

# ${var : start-here : go-this-many}
#echo ${res:0:1}
#echo ${res:1:1}

[ ${res:0:1} == '-' ] && _err "must be positive number"
[ ${res:0:1} == "0" ] && [ ${res:1:1} != '.' ] && _err "no leading 0's"
[ ${res:0:1} == "0" ] && [ ${res:1:2} == '.0' ] && _err "number is < 0.1"

i=0;dots=0;digits=0
echo && while [ "${res:$i:1}" != '' ]
do
    printf " %s" ${res:$i:1}
    [ "${res:$i:1}" == \. ] && ((dots++))
    [ $dots -eq 0 ] && ((digits++))
    ((i++))
    [ $i -gt 5 ] && _err 'number is too long'
done && echo
[ $digits -gt 2 ] && _err 'number is >= 100'
[ $dots -gt 1 ] && _err 'too many decimal points'


echo -e "\n Setting resolution to  -->  $res % \n"

path=`pwd`
dirs="MCbetas MCgammas"

file1="proba.car"

for dir in $dirs
do

    here="$path/$dir"

    sed -i -c /scriptme/s/res\ =\ .*\ \!scriptme/res\ =\ ${res}\ \!scriptme/ "$here/$file1"

done

exit

