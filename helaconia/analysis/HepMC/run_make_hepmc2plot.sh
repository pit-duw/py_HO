#!/usr/bin/env bash

echo "Please enter in the plot file with .f90 extension :"
while read FORFILE
do
    if [ ! -f "$FORFILE" ]
    then
        echo "$FORFILE does not exist"
        echo "exit(y/n) ?"
        read x
        if [ "$x" == "y" ]
        then
            exit
        else
            echo "input the valid plot file with .f90 extension again"
        fi
    else
        break
    fi
done
echo "The HepMC plot file is $FORFILE"
arr=$(echo $FORFILE | tr ".f90" "\n")
FORFILE=$arr.o
make clean -f makefile_hepmc2plot
make -f makefile_hepmc2plot FORFILE=$FORFILE LIBGFORTRANPATH="-L/usr/local/gfortran/lib" LIBGFORTRAN=-lgfortran EXTRALIBS="-lfjcore -lheptoptagger" EXTRAPATHS="-L../../lib" INCLOPTION="" HEPMCLOCATION="/Users/erdissshaw/Works/HepMC/install"
 
