#!/bin/bash

echo "Please enter in the C++ HepMC2Root main file :"
while read MAINFILE
do
    if [ ! -f "$MAINFILE" ]
    then
        echo "$MAINFILE doesnot exist"
        echo "exit(y/n) ?"
        read x
        if [ "$x" == "y" ]
        then
            exit
        else
            echo "input the valid C++ HepMC2Root main file again"
        fi
    else
        break
    fi
done
echo "The C++ HepMC2Root main file is $MAINFILE"

make clean -f makefile_hepmc2root
make -f makefile_hepmc2root MAINFILE=$MAINFILE FORFILE=$FORFILE LIBGFORTRANPATH= LIBGFORTRAN= EXTRALIBS="-lfjcore -lheptoptagger" EXTRAPATHS= INCLOPTION="-I../../jets/fjcore -I../heptoptagger" HEPMCLOCATION=/Users/erdissshaw/Works/HepMC/install ROOTLOCATION=/Users/erdissshaw/Works/ROOT/root
# ROOT LIBS
# -lTree -lCore -lGpad -lHist -lGraf -lGraf3d -lRint -lPostscript -lMatrix -lMathCore -lRIO -lNet -lThread -lCint -lEG
