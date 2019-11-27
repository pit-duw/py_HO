#!/usr/bin/env bash

echo "Please enter in the analysis file with extension .f :"
while read PROCESS
do
    if [ ! -f "$PROCESS" ]
    then
        echo "$PROCESS doesnot exist"
        echo "exit(y/n) ?"
        read x
        if [ "$x" == "y" ]
        then
            exit
        else
            echo "input the valid analysis file with extension .f again"
        fi
    else
        break
    fi
done
echo "The analysis file is $PROCESS"
arr=$(echo $PROCESS | tr ".f" "\n")
PROCESS=$arr.o
make clean -f makefile_lhe2topdrawer
make -f makefile_lhe2topdrawer PROCESS=$PROCESS

