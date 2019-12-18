#!/bin/bash
# wrapper for the functions which run PYTHIA8
thisdir=`pwd`

function runShower {

    compileShower

    UMASSPY8=$UMASS
    DMASSPY8=$DMASS
    SMASSPY8=$SMASS
    CMASSPY8=$CMASS
    if (( $(bc <<< "$BMASS < 0.0") )); then BMASSPY8=$BMASS; fi
    if (( $(bc <<< "$BMASS >= 0.0") )); then BMASSPY8=$BMASS; fi
    EMASSPY8=$EMASS
    MUMASSPY8=$MUMASS
    TAUMASSPY8=$TAUMASS
    GMASSPY8=$GMASS

    if [ $BEAM1 == 1 ]; then 
	iA=2212; 
    elif [ $BEAM1 == -1 ]; then
	iA=-2212;
    elif [ $BEAM1 == 0 ]; then
	iA=-11;
    elif [ $BEAM1 == 2 ]; then
	iA=2112;
    elif [ $BEAM1 == -2 ]; then
	iA=-2112;
    else
	echo "ERROR:Unknown incoming beam 1";
	exit 1;
    fi
    
    if [ $BEAM2 == 1 ]; then
	iB=2212;
    elif [ $BEAM2 == -1 ]; then
	iB=-2212;
    elif [ $BEAM2 == 0 ]; then
	iB=11;
    elif [ $BEAM2 == 2 ]; then
	iB=2112;
    elif [ $BEAM2 == -2 ]; then
	iB=-2112;
    else
	echo "ERROR:Unknown incoming beam 2";
	exit 1;
    fi

    whichpytpdf $PYTPDF
    if [ $pdftype == 1 ] # PYTPDF == EXTPDF or extpdf
    then
	whichpdflib $PDFLIBRARY
	if [ $UsedPdfLib == "THISLIB" ]
	then
	    dothelinks
	elif [ $UsedPdfLib == "LHAPDF" ]
	then
	    linklhapdf
	fi
    fi
    NEVENTS_TOT=$(grep -rin "<event>" $EVTFILE | wc -l)
    if [ "$NEVENTS"  -lt 0 ] || [ "$NEVENTS" -gt "$NEVENTS_TOT" ];
    then
	NEVENTS=$NEVENTS_TOT
    fi
    ERR_NUM_PY8=$(echo "$ERR_NUM_PY8 * $NEVENTS" | bc)
    ERR_NUM_PY8=${ERR_NUM_PY8/.*}
    let "ERR_NUM_PY8 += 1"
    echo " "
    echo "***** Now ready for showering" $NEVENTS "events with PYTHIA8 "
    echo " "

    runHOPY8
    teststatus runShower
}

function compileShower {

    if [ "$PDFCODE" -ne 0 ]
    then
	PYTPDF='EXTPDF'
	PDFLIBRARY='LHAPDF'
	UsedPdfLib='LHAPDF'
	export LIBLHAPDF="-L$LHAPDFPATH/lib -lLHAPDF"
	LHALINK=DYNAMIC
	if [ "$PDFCODE" -gt 0 ]; then LHAOFL=FREEZE; fi
	if [ "$PDFCODE" -gt 1 ]; then PDFSET=$PDFCODE; fi
	if [ "$PDFCODE" -lt 0 ]; then LHAOFL=EXTRAPOLATE; fi
	if [ "$PDFCODE" -lt -1 ]; then PDFSET=$((PDFCODE*-1));fi
	linklhapdf
	# convert PDF number to name reading PDFsets.index
	while read pdfline
	do
	    pdfinfo=($(echo ${pdfline}))
	    if [ "$PDFSET" == "${pdfinfo[0]}" ]
	    then
		PDFSETN="${pdfinfo[4]}"
	    fi
	done < $LHAPDFPATH/share/lhapdf/PDFsets.index
    elif [ "$PDFCODE" -eq 0 ]
    then
	PYTPDF='DEFAULT'
	PDFLIBRARY='THISLIB'
	UsedPdfLib='THISLIB'
	export LIBLHAPDF=-llhapdfdummy # it is in Pythia8 lib
	# the following is dummy
	LHALINK=DYNAMIC
	LHAOFL=FREEZE
    else
	echo 'Invalid PDFCODE' $PDFCODE
	exit 1
    fi
    
    export PYTHIA8LOCATION=$PY8PATH

    if [ -f $thisdir/Pythia8.exe ]
    then
	rm -rf $thisdir/Pythia8.exe
    fi

    if [ -d $thisdir/temp ]
    then
	rm -rf $thisdir/temp
    fi
    mkdir $thisdir/temp

    nfile=0
    nf90file=0
    ncppfile=0
    if [ "$PY8UTI" == "" ]
    then
	cp -rf $thisdir/PYTHIA8/Pythia8_lhe2hep.cc $thisdir/temp/Pythia8.cc
	cp -rf $thisdir/PYTHIA8/makefile_lhe2hep $thisdir/temp/makefile
    else
	if [ "$MLM_MERGING" == "1" ]
	then
	    path_insert="/MLMmerging"
	else
	    path_insert=""
	fi
	arr=$(echo $PY8UTI | tr " " "\n")
	#nfile=0
	#nf90file=0
	#ncppfile=0
	for xx in $arr
	do
	    if [ -f "$HODirectory/analysis$path_insert/PYTHIA8/$xx.f90" ]
	    then
		nfile=`expr $nfile + 1`
		nf90file=`expr $nf90file + 1`
		PY8UTIF90[$nfile]="$xx.f90"
		PY8UTIOBJ[$nfile]="$xx.o"
	    elif [ -f "$HODirectory/analysis$path_insert/PYTHIA8/$xx.cc" ]
	    then
		nfile=`expr $nfile + 1`
		ncppfile=`expr $ncppfile + 1`
		PY8UTICPP[$nfile]="$xx.cc"
		PY8UTIHF[$nfile]="$xx.h"
		PY8UTIOBJ[$nfile]="$xx.o"
	    else
		echo "WARNING:Cannot find $xx.f90 or $xx.cc! Ignore it!"
	    fi
	done
	nf90cpp=`expr $ncppfile \* $nf90file`
	if [ $nfile -gt 0 ] && [ "$nf90cpp" == "0" ] 
	then
	    if [ "$MLM_MERGING" == "0" ]
	    then
		if [ $nf90file -gt 0 ]
		then
		    cp -rf $thisdir/PYTHIA8/Pythia8_lhe2plot.cc $thisdir/temp/Pythia8.cc
		    cp -rf $thisdir/PYTHIA8/makefile_lhe2plot $thisdir/temp/makefile
		    cp -rf $thisdir/various/LHEFRead.h $thisdir/temp/LHEFRead.h
		    cp -rf $thisdir/various/LHEF.h $thisdir/temp/LHEF.h
		else
		    cp -rf $thisdir/PYTHIA8/makefile_lhe2root $thisdir/temp/makefile
		fi
	    else
                if [ $nf90file -gt 0 ]
                then
                    cp -rf $thisdir/PYTHIA8/Pythia8_MLM_lhe2plot.cc $thisdir/temp/Pythia8.cc
                    cp -rf $thisdir/PYTHIA8/makefile_MLM_lhe2plot $thisdir/temp/makefile
                    cp -rf $thisdir/various/LHEFRead.h $thisdir/temp/LHEFRead.h
                    cp -rf $thisdir/various/LHEF.h $thisdir/temp/LHEF.h
		    cp -rf $thisdir/various/MLM_merge.f90 $thisdir/temp/MLM_merge.f90
		    cp -rf $thisdir/various/HEPEUPf90.inc $thisdir/temp/HEPEUPf90.inc
                #else
                #    cp -rf $thisdir/PYTHIA8/makefile_lhe2root $thisdir/temp/makefile
                fi
	    fi
	    # reset PY8UTI
	    PY8UTI=""
	    for ii in $(seq 1 $nfile)
	    do
		if [ $nf90file -gt 0 ]
		then
		    cp -rf $HODirectory/analysis$path_insert/PYTHIA8/${PY8UTIF90[$ii]} $thisdir/temp/${PY8UTIF90[$ii]}
		    if [ $ii == 1 ]
                    then
			PY8UTI="${PY8UTIOBJ[$ii]}"
                    else
			PY8UTI="$PY8UTI ${PY8UTIOBJ[$ii]}"
                    fi
		    if [ "$MLM_MERGING" == "1" ]
		    then
			PY8UTI="$PY8UTI MLM_merge.o"
		    fi
		else
		    cp -rf $HODirectory/analysis$path_insert/PYTHIA8/${PY8UTICPP[$ii]} $thisdir/temp/Pythia8.cc
		    if [ -f "$HODirectory/analysis/PYTHIA8/${PY8UTIHF[$ii]}" ]
		    then
			cp -rf $HODirectory/analysis$path_insert/PYTHIA8/${PY8UTIHF[$ii]} $thisdir/temp/${PY8UTIHF[$ii]}
		    fi
		fi
		#if [ $ii == 1 ]
		#then
		#    PY8UTI="${PY8UTIOBJ[$ii]}"
		#else
		#    PY8UTI="$PY8UTI ${PY8UTIOBJ[$ii]}"
		#fi
	    done
	    if [ $nf90file -gt 0 ]
	    then
		cp -rf $HODirectory/analysis/include/HEPMC90.INC $thisdir/temp/HEPMC90.INC
		cp -rf $HODirectory/analysis/hbook/dbookf90.inc $thisdir/temp/dbookf90.inc
	    fi
	else
	    if [ $ncppf90 -gt 0 ]
	    then
		echo "WARNING:Two types, i.e. lhe2plot and lhe2root"
		echo "WARNING:Take the lhe2hep mode"
	    fi
	    cp -rf $thisdir/PYTHIA8/Pythia8_lhe2hep.cc $thisdir/temp/Pythia8.cc
	    cp -rf $thisdir/PYTHIA8/makefile_lhe2hep $thisdir/temp/makefile
	    # reset PY8UTI
	    PY8UTI=""
	    nfile=0
	    ncppfile=0
	    nf90file=0
	fi
    fi
    
    cp -rf $PYTHIA8LOCATION/xmldoc $thisdir/temp
    cp -rf $PYTHIA8LOCATION/xmldoc $thisdir
    cd $thisdir/temp
    PY8UTIMAKE="PY8UTI="$PY8UTI
    FORFILEMAKE="FORFILE="$PY8UTI
    if [ "$PY8UTI" == "" ] && [ "$nfile" == "0" ]
    then
        LGFORTRANMAKE="LIBGFORTRAN=-lgfortran"
        LGFORTRANPATHMAKE="LIBGFORTRANPATH="$LIBGFORTRANPATH
    else
	if [ $nf90file -gt 0 ]
	then
	# try to find the correct path for libgfortran
	# from PATH, LD_LIBRARY_PATH, DYLD_LIBRARY_PATH
	    arr=$(echo $PATH | tr ":" "\n")
	    LIBGFORTRANPATH=""
	    for xx in $arr
	    do
		if [ -f "$xx/libgfortran.a" ]
		then
		    LIBGFORTRANPATH="-L$xx"
		    break
		fi
	    done
	    if [ "$LIBGFORTRANPATH" == "" ]
	    then
		arr=$(echo $LD_LIBRARY_PATH | tr ":" "\n")
		for xx in $arr
		do
		    if [ -f "$xx/libgfortran.a" ]
		    then
			LIBGFORTRANPATH="-L$xx"
			break
		    fi
		done
		if [ "$LIBGFORTRANPATH" == "" ]
		then
		    arr=$(echo $DYLD_LIBRARY_PATH | tr ":" "\n")
		    for xx in $arr
		    do
			if [ -f "$xx/libgfortran.a" ]
			then
			    LIBGFORTRANPATH="-L$xx"
			    break
			fi
		    done
		fi
	    fi
	    if [ "$LIBGFORTRANPATH" == "" ]
	    then
		echo "WARNING:Cannot find gfortran lib in PATH,LD_LIBRARY_PATH and DYLD_LIBRARY_PATH"
	    fi
	    LGFORTRANMAKE="LIBGFORTRAN=-lgfortran"
            LGFORTRANPATHMAKE="LIBGFORTRANPATH="$LIBGFORTRANPATH
	else
	    LGFORTRANMAKE="LIBGFORTRAN="
	    LGFORTRANPATHMAKE="LIBGFORTRANPATH="
	    if [ -r "$HODirectory/input/paths/rootpath" ]
	    then
		read -r rootpath<"$HODirectory/input/paths/rootpath"
	    fi
	fi
	if [ -r "$HODirectory/input/paths/fastjetpath" ]
	then
	# link to fastjet
	    read -r fastjetpath<"$HODirectory/input/paths/fastjetpath"
	    fastjetlibflag=`$fastjetpath --libs`
	    ExtraLibs="$ExtraLibs $fastjetlibflag"
	    if [[ $ExtraLibs != *"-lfastjetplugins"* ]]
	    then
		ExtraLibs="$ExtraLibs -lfastjetplugins"
	    fi
	    if [[ $ExtraLibs != *"-lsiscone"* ]]
	    then
		ExtraLibs="$ExtraLibs -lsiscone -lsiscone_spherical"
	    fi
	    if [[ $ExtraLibs != *"-lfastjetfortran"* ]]
	    then
		ExtraLibs="$ExtraLibs -lfastjetfortran"
	    fi
	    if [[ $ExtraPaths != *"-L\$(HOPATH)/lib"* ]]
	    then
		ExtraPaths="$ExtraPaths -L\$(HOPATH)/lib"
	    fi
	    if [[ $ExtraLibs != *"-lheptoptagger"* ]]
	    then
		if [ -f "$HODirectory/lib/libheptoptagger.a" ]
		then
		    ExtraLibs="$ExtraLibs -lheptoptagger"
		else
		    echo "WARNING:Failed to compile heptoptagger"
		    echo "WARNING:The functionality will be switched off"
		    echo "WARNING:Possible reasons are:"
		    echo "WARNING:1) Cannot find libgfortran"
		    echo "WARNING:2) the version of gcc >= 4.5"
		fi
	    fi
	    if [ $ncppfile -gt 0 ]
	    then
		fastjetbinpath=`dirname $fastjetpath`
		fastjetmainpath=`dirname $fastjetbinpath`
		fastjetlibpath=$fastjetmainpath/lib
		fastjetincludepath=$fastjetmainpath/include
		ExtraIncl="$ExtraIncl -I$fastjetincludepath"
		if [ -f "$HODirectory/lib/libheptoptagger.a" ]
		then
		    ExtraIncl="$ExtraIncl -I\$(HOPATH)/analysis/heptoptagger"
		fi
	    fi
	else
	# use fjcore, which cannot use jet area and siscone
	    if [[ $ExtraLibs != *"-lfjcore"* ]]
	    then
		ExtraLibs="$ExtraLibs -lfjcore"
	    fi
	    if [[ $ExtraPaths != *"-L\$(HOPATH)/lib"* ]]
	    then
		ExtraPaths="$ExtraPaths -L\$(HOPATH)/lib"
	    fi
	    if [[ $ExtraLibs != *"-lheptoptagger"* ]]
	    then
                if [ -f "$HODirectory/lib/libheptoptagger.a" ]
                then
                    ExtraLibs="$ExtraLibs -lheptoptagger"
                else
                    echo "WARNING:Failed to compile heptoptagger"
                    echo "WARNING:The functionality will be switched off"
                    echo "WARNING:Possible reasons are:"
                    echo "WARNING:1) Cannot find libgfortran"
                    echo "WARNING:2) the version of gcc >= 4.5"
                fi
	    fi
	    if [ $ncppfile -gt 0 ]
            then
                ExtraIncl="$ExtraIncl -I\$(HOPATH)/jets/fjcore"
                if [ -f "$HODirectory/lib/libheptoptagger.a" ]
		then
                    ExtraIncl="$ExtraIncl -I\$(HOPATH)/analysis/heptoptagger"
                fi
            fi
	fi
    fi
    # include the path of HepMC if it exsits
    if [ "$HepMCDirectory" != "" ] && [ -d "$HepMCDirectory/lib" ] && [ "$ncppfile" == "0" ]
    then
	ExtraPaths="$ExtraPaths -L$HepMCDirectory/lib"
	echo "INFO:If you encounter the problem of finding library libHepMC"
	echo "INFO:Please also add HepMC to LD_LIBRARY_PATH/DYLD_LIBRARY_PATH and PATH in .bashrc"
    fi
    LIBSMAKE="EXTRALIBS="$ExtraLibs
    LIBSMAKEP="EXTRAPATHS="$ExtraPaths
    INCLMAKE="INCLOPTION="$ExtraIncl
    HODIRMAKE="HOPATH="$HODirectory
    ROOTLOCATIONMAKE="ROOTLOCATION="$rootpath
    echo Pythia8 "$PY8UTIMAKE" "$FORFILEMAKE" "$LGFORTRANPATHMAKE" "$LGFORTRANMAKE" "$LIBSMAKE" "$LIBSMAKEP" "$INCLMAKE" "$HODIRMAKE" "$ROOTLOCATIONMAKE"
    make Pythia8 "$PY8UTIMAKE" "$FORFILEMAKE" "$LGFORTRANPATHMAKE" "$LGFORTRANMAKE" "$LIBSMAKE" "$LIBSMAKEP" "$INCLMAKE" "$HODIRMAKE" "$ROOTLOCATIONMAKE"

    cd $thisdir
    if [ ! -f $thisdir/temp/Pythia8.exe ]
    then
	echo "Pythia8 compilation did not succeed, exiting"
	exit -1
    fi
    cp $thisdir/temp/Pythia8.exe $thisdir
    #rm -rf $thisdir/temp

    teststatus compileShower
}

function runHOPY8 {
    ifile="HELACOnia_PYTHIA8_input"
    #exefile="HELACOnia_PYTHIA8_EXE"
    if [ -f ./$ifile ]
    then
	rm ./$ifile
    fi
    cat <<EOF > $ifile

! 1) Settings used in the main program.
Main:numberOfEvents = $NEVENTS_TOT     ! Number of events in the LHE file
Main:spareMode1 = $NEVENTS             ! Number of events to be showered
Main:spareWord1 = $EVENT_NORM          ! Event weights are normalized to sum
                                       ! or average to the cross section
Main:timesAllowErrors = $ERR_NUM_PY8   ! Number of allowed errors
Main:showChangedSettings = on          ! Shows all non-default settings
Main:showChangedParticleData = off     ! Shows all non-default particle settings

! 2) Settings related to output in init(), next(), and stat().
Init:showChangedSettings = on          ! Shows all non-default settings
Init:showChangedParticleData = off     ! Shows all non-defulat particle settings
Next:numberCount = 100                 ! print message every n events
Next:numberShowInfo = 2                ! print event information n times
Next:numberShowProcess = 1             ! print process record n times
Next:numberShowEvent = $MAXPR_PY8      ! print event record n times
Stat:showProcessLevel = on             ! Process statistics
Stat:showErrors = on                   ! Error statistics
Check:epTolErr = 0.001                 ! Momentum-conservation tolerance

! 3) Beam-parameter settings.
Beams:idA = $iA                        ! Beam identities
Beams:idB = $iB                        ! Beam identities
Beams:frameType = 4                    ! LHE initialization
Beams:LHEF = $EVTFILE               ! Input LHE file

! 4) Switch on/off the key event-generation steps.
EOF

    if [ $UsedPdfLib = "LHAPDF" ]
    then
	cat <<EOF >> $ifile
PDF:useLHAPDF = on                     ! Use of LHAPDF
PDF:LHAPDFSET = $PDFSETN               ! PDF set
EOF
	if [ $PDFGROUP = "LHAPDF" ]
	then
	    cat <<EOF >> $ifile
PDF:extrapolateLHAPDF = off            ! extrapolate PDF set outside the boundaries
EOF
	elif [ $PDFGROUP = "LHAEXT" ]
	then
	    cat <<EOF >> $ifile
PDF:extrapolateLHAPDF = on             ! extrapolate PDF set outside the boundaries
EOF
	else
	    echo "ERROR:Unknown PDFGROUP " $PDFGROUP
	    exit 1
	fi
    else
	cat <<EOF >> $ifile
PDF:useLHAPDF = off                    ! Use of LHAPDF
PDF:pSet = 7                           ! CTEQ6L
EOF
    fi

    if ((`bc <<< "$LAMBDAPYTH >= 0.0"`))
    then
	cat <<EOF >> $ifile
#coupSM:Lambda5 = $LAMBDAPYTH          ! five-flavour lLambdaQCD
EOF
    fi
    cat <<EOF >> $ifile
ProcessLevel:all = on                  ! Generation
ProcessLevel:resonanceDecays = on      ! Resonance decays
PartonLevel:all = on                   ! Parton level: if off, stops after hard process generation
EOF
    if [ $UE_PY8 = ".FALSE." ]
    then
	cat <<EOF >> $ifile
PartonLevel:MPI = off                  ! Multiple interactions
EOF
    fi
    
    cat <<EOF >> $ifile
PartonLevel:ISR = on                   ! Initial state shower
PartonLevel:FSR = on                   ! Final state shower
PartonLevel:FSRinProcess = on          ! Final state shower in association with the hard process
PartonLevel:FSRinResonances = on       ! Final state shower in resonance decays
HadronLevel:all = on                   ! Hadron level: if off, stops before hadronization
EOF
    
    if [ $HADRONIZE_PY8 = ".FALSE." ]
    then
	cat <<EOF >> $ifile
HadronLevel:Hadronize = off            ! Hadronization
EOF
    else
	cat <<EOF >> $ifile
HadronLevel:Hadronize = on             ! Hadronization
EOF
    fi

    cat <<EOF >> $ifile
#HadronLevel:Decay = on                ! Hadron decays
PhaseSpace:mHatMin = 4.                ! Min invariant mass
PhaseSpace:mHatMax = -1.               ! Max invariant mass
PhaseSpace:pTHatMin = 0.               ! Min pT in 2->2
PhaseSpace:pTHatMax = -1.              ! Max pT in 2->2
PhaseSpace:pTHatMinDiverge = 1.        ! If massless final state, to avoid divergences
PhaseSpace:useBreitWigners = on        ! Masses according to Breit-Wigner
#PhaseSpace:pTHat3Min = 0.             ! Min pT for the hardest parton in 2->3
PhaseSpace:pTHat3Max = -1.             ! Max pT for the hardest parton in 2->3
PhaseSpace:pTHat5Min = 0.              ! Min pT for the softest parton in 2->3
PhaseSpace:pTHat5Max = -1.             ! Max pT for the softest parton in 2->3
PhaseSpace:RsepMin = 0.                ! Min R separation in 2->3

! 5) Final-state shower.
TimeShower:pTmaxMatch = 1              ! Use scalup (re-check)
TimeShower:pTmaxFudge = 1.             ! Factor changing the max scale
TimeShower:alphaSvalue = 0.118         ! Alpha_s(MZ) in final-state shower
TimeShower:alphaSorder = 1             ! Alpha_s running order in final-state shower
TimeShower:alphaEMorder = 0            ! Alpha_EM running order in final-state shower
TimeShower:interleave = on             ! If on, FSR interleaved with ISR
TimeShower:allowBeamRecoil = on        ! If off, no energy transfer from ISR to FSR
TimeShower:dampenBeamRecoil = off      ! Dampens the effect of beam recoil
TimeShower:globalRecoil = on           ! All final-state particles recoil against the branching
TimeShower:nMaxGlobalRecoil =  1       ! Number of splittings with TimeShower:globalRecoil = on
TimeShower:globalRecoilMode = 2        ! Global recoil only for S events whose first emission is FSR
TimeShower:nMaxGlobalBranch = 1        ! Number of FSR splittings proposed with global recoil
TimeShower:nPartonsInBorn = -1         ! Number of Born QCD final-state partons (to treat H and S differently)
TimeShower:limitPTmaxGlobal = on       ! Limits pT < min(SCALUP,mDipole/2)
TimeShower:QCDshower = on              ! QCD final-state shower
TimeShower:nGluonToQuark = 5           ! Number if flavors allowed in g->qqbar
TimeShower:QEDshowerByQ = off          ! Prevent quarks from radiating photons
TimeShower:QEDshowerByL = off          ! Prevent leptons from radiating photons
TimeShower:QEDshowerByGamma = off      ! Prevent photons from branching
TimeShower:MEcorrections = off         ! No Matrix-element corrections
TimeShower:MEafterFirst = off          ! No Matrix-element corrections after first emission
TimeShower:phiPolAsym = on             ! Azimuthal asymmetry induced by gluon polarization
TimeShower:alphaSuseCMW = false        ! Use the CMW prescription in FSR

! 6) Initial-state shower.
SpaceShower:pTmaxMatch = 1             ! Use scalup (re-check)
SpaceShower:pTmaxFudge = 1.            ! Factor changing the max scale
SpaceShower:alphaSvalue = 0.118        ! Alpha_s(MZ) in initial-state shower
SpaceShower:alphaSorder = 1            ! Alpha_s running order in initial-state shower
SpaceShower:alphaEMorder = 0           ! Alpha_EM running order in initial-state shower
SpaceShower:QCDshower = on             ! QCD initial-state shower
SpaceShower:QEDshowerByQ = off         ! Prevent quarks from radiating photons
SpaceShower:QEDshowerByL = off         ! Prevent leptons from radiating photons
SpaceShower:MEcorrections = off        ! No Matrix-element corrections
SpaceShower:MEafterFirst = off         ! No Matrix-element corrections after first emiision
SpaceShower:phiPolAsym = on            ! Azimuthal asymmetry induced by gluon polarization
SpaceShower:nQuarkIn = 5               ! Number of flavors in g->qqbar and also in incoming beams
SpaceShower:rapidityorder = off        ! Do not order branchings in rapidity
SpaceShower:alphaSuseCMW = false       ! Use the CMW prescription in ISR

! 7) Non-perturbative stuff
BeamRemnants:primordialKT = off        ! No primordial kT

! 8) Particle characteristics.
1:m0 = $DMASSPY8                       ! down mass
2:m0 = $UMASSPY8                       ! up mass
3:m0 = $SMASSPY8                       ! strange mass
4:m0 = $CMASSPY8                       ! charm mass
5:m0 = $BMASSPY8                       ! bottom mass
6:m0 = $TMASS                          ! top mass
11:m0 = $EMASSPY8                      ! electron mass
13:m0 = $MUMASSPY8                     ! muon mass
15:m0 = $TAUMASSPY8                    ! tauon mass
23:m0 = $ZMASS                         ! Z mass
24:m0 = $WMASS                         ! W mass
25:m0 = $HGGMASS                       ! Higgs mass
6:mWidth = $TWIDTH                     ! top width
23:mWidth = $ZWIDTH                    ! Z width
24:mWidth = $WWIDTH                    ! W width
25:mWidth = $HGGWIDTH                  ! Higgs width 
EOF
    
    if [ $PI_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
111:mayDecay = false                   ! stable pi0
#211:mayDecay = false                  ! stable pions
EOF
    fi
    if [ $B_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
511:maydecay = false                   ! stable B hadrons
521:maydecay = false                   ! stable B hadrons
531:maydecay = false                   ! stable B hadrons
541:maydecay = false                   ! stable B hadrons
553:maydecay = false                   ! stable B hadrons
5112:maydecay = false                  ! stable B hadrons
5122:maydecay = false                  ! stable B hadrons
5132:maydecay = false                  ! stable B hadrons
5222:maydecay = false                  ! stable B hadrons
5232:maydecay = false                  ! stable B hadrons
5332:maydecay = false                  ! stable B hadrons
EOF
    fi

    if [[ $WP_STABLE_PY8 = ".TRUE." || $WM_STABLE_PY8 = ".TRUE." ]]
    then
	cat <<EOF >> $ifile
24:maydecay = false                    ! stable W boson
EOF
    fi
    if [ $Z_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
23:maydecay =  false                   ! stable Z boson
EOF
    fi
    if [ $H_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
25:maydecay = false                    ! stable Higgs boson
EOF
    fi
    if [ $TAUP_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
-15:maydecay = false                   ! stable tau+
EOF
    fi
    if [ $TAUM_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
15:maydecay = false                    ! stable tau-
EOF
    fi
    if [ $MUP_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
-13:maydecay = false                   ! stable mu+
EOF
    fi
    if [ $MUM_STABLE_PY8 = ".TRUE." ]
    then
	cat <<EOF >> $ifile
13:maydecay = false                    ! stable mu-
EOF
    fi
    if [ $RNDEVSEED_PY8 != 0 ]
    then
	cat <<EOF >> $ifile
Random:setSeed = on                    ! random seed
Random:seed = $RNDEVSEED_PY8
EOF
    fi
    # following is for the decay chains
    #arr=()
    #numDM=0
    #iDM=0
    #for i in {1..99}
    #do
    #eval arr_elem='$'DM_$i
#	if [ "$arr_elem" !="" ]
#	then
#	    arr=("${arr[@]}" "$arr_elem")
#	fi
#    done
    
#    numDM=${#arr[@]}
#    for DM in "${arr[@]}"
#    do
#	DM=( $DM )
#	cat <<EOF >> ./$ifile
#${DM[0]}:onMode = off                     ! First switch decays off
#EOF
#    done

#    for DM in "${arr[@]}"
#    do
#	let "iDM+=1"
#	DM=( $DM )
#	num_of_part=$(echo ${#DM[@]})
#	let "num_of_part-=6"

#	if [ $num_of_part = 2 ]
#	then
#	    cat <<EOF >> ./$ifile
#${DM[0]}:onIfMatch = ${DM[2]} ${DM[3]}      ! Decay mode $iDM of $numDM
#EOF
#	elif [ $num_of_part = 3 ]
#	then
#	    cat <<EOF >> ./$ifile
#${DM[0]}:onifMatch = ${DM[2]} ${DM[3]} ${DM[4]}  ! Decay mode $iDM of $numDM
#EOF
#	elif [ $num_of_part = 4 ]
#	then
#	    cat <<EOF >> ./$ifile
#${DM[0]}:onifMatch = ${DM[2]} ${DM[3]} ${DM[4]} ${DM[5]} ! Decay mode $iDM of $numDM
#EOF
#	elif [ $num_of_part = 5 ]
#	then
#	    cat <<EOF >> ./$ifile
#${DM[0]}:onifMatch = ${DM[2]} ${DM[3]} ${DM[4]} ${DM[5]} ${DM[6]} ! Decay mode $iDM of $numDM
#EOF
#	fi
#    done
    cp -f $thisdir/$ifile Pythia8_lhe.cmnd
    rm $thisdir/$ifile
}
	

# this function set the parameter pdftype according to the value
# of PYTPDF (the entry of this function) given in input
function whichpytpdf {
    case $1 in
	DEFAULT|default) pdftype=0 ;;
	EXTPDF|extpdf) pdftype=1 ;;
	*) echo "ERROR:no such option in whichpytpdf"; exit 1;;
    esac
}

# checks that the value given to PDFLIBRARY in input is meaningful
function whichpdflib {
    case $1 in
	THISLIB|thislib) UsedPdfLib=THISLIB ;;
	PDFLIB|pdflib) UsedPdfLib=PDFLIB ;;
	LHAPDF|lhapdf) UsedPdfLib=LHAPDF ;;
	*) echo "ERROR:no such library for PDFS; failure in whichpdflib"; exit 1 ;;
    esac
}

# checks that the value given to LHALINK in input is meaningful
function whichlhapdf {
    case $1 in 
	STATIC|static) UsedLhaPdf=lhasta ;;
	DYNAMIC|dynamic) UsedLhaPdf=lhadyn ;;
	*) echo "ERROR:no such option for LHAPDF; failure in whichlhapdf"; exit 1;;
    esac
}

# $? is the value of last exectued command. A call to this function
# after a failure will cause the program to quit the script
function teststatus {
    rc=$?
    if [ 0 = $rc ]
    then
	:
    else
	echo $* did not succed, exit status=$rc
	exit $rc
    fi
}

# utility function for dothelinks
# remove extension for example haha.txt -> haha
function stripextension {
    echo $1 | sed "s/\..*\$//"
}

# utility function for dothelinks
# to be capital
function capitalize {
    echo $1 | sed "y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/"
}

# creates logical links for the PDF grid files. By P. Nason
function dothelinks {
    if [ -d $PDFPATH ]
    then
	for i in ${PDFPATH}/*.dat ${PDFPATH}/*.tbl
	do
	    if [ -f $i ]
	    then
		name=`basename $i`
		name=`stripextension $name`
		case $name in
		    mrst200*) ;;
		    *mrs*|grpol|grvdm|lac1|pion[1-3]) name=`capitalize $name`;;
		esac
		if [ ! -L $thisdir/$name ] || [ $thisdir/$name -ef $i ]
		then
		    ln -sf $i $thisdir/$name
		fi
	    fi
	done
	for i in ${PDFPATH}/a02.*
	do
	    if [ -f $i ]
	    then
		name=`basename $i`
		if [ ! -L $thisdir/$name ] || [ ! $thisdir/$name -ef $i ]
		then
		    ln -sf $i $thisdir/$name
		fi
	    fi
	done
    fi
}

# creates logical links for LHAPDF, and replaced PDF group name (unused
# by LHAPDF) with a LHAPDF-specific string
function linklhapdf {
    case $LHAOFL in
	FREEZE|freeze) PDFGROUP=LHAPDF ;;
	EXTRAPOLATE|extrapolate) PDFGROUP=LHAEXT ;;
	*) echo "ERROR:no such option; failure in linklhapdf"; exit 1 ;;
    esac
    #source ../Source/fj_lhapdf_opts
}
		    