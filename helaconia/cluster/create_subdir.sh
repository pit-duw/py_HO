#!/bin/bash
if [ $# -ge 1 ]
then
    subdir=$1
    if [ $# -ge 2 ]
    then
	addonfile=$2
    fi
    #force = $2
else
    echo "1 argument is necessary"
    exit
fi

cd ..
hodir=$PWD
cd cluster
# global flags
flag_dir="n"
flag_dir_x="n"
flag_cpfile="n"
flag_cpfile_x="n"
flag_lnfile="n"
flag_lnfile_x="n"
# helper functions
create_dir(){
    local dirname="$1"
    if [ -d "$dirname" ]
    then
	if [ "$flag_dir" == "y" ]
	then
	    if [ "$flag_dir_x" == "y" ]
	    then
		rm -rf $dirname
		mkdir $dirname
	    else
		return
	    fi
	fi
	echo "$dirname already exists"
	echo "do you want to create a new one ? (y/n)"
	read local x
	echo "do you want to force this option to all directories (y/n)"
	read flag_dir
	if [ "$flag_dir" == "y" ]
	then
	    flag_dir_x=$x
	fi
	if [ "$x" == "y" ]
	then
	    rm -rf $dirname
	    mkdir $dirname
	fi
    else
	mkdir $dirname
    fi
}
is_file_exists(){
    local fname="$1"
    if [ -f "$fname" ]
    then
	return 1
    else
	return 0
    fi
}
copy_file(){
    local newfname="$2"
    local oldfname="$1"
    if [ ! -f "$oldfname" ]
    then
	echo "ERROR:Cannot find $oldfname"
	exit 1
    fi
    if [ -f "$newfname" ]
    then
	if [ "$flag_cpfile" == "y" ]
        then
            if [ "$flag_cpfile_x" == "y" ]
	    then
		rm $newfname
		cp $oldfname $newfname
            else
                return
            fi
        fi
	echo "$newfname already exists"
        echo "do you want to replace it ? (y/n)"
        read local x
	echo "do you want to force this option to all copied files (y/n)"
        read flag_cpfile
        if [ "$flag_cpfile" == "y" ]
        then
            flag_cpfile_x=$x
        fi
        if [ "$x" == "y" ]
        then
            rm $newfname
            cp $oldfname $newfname
        fi
    else
	cp $oldfname $newfname
    fi
}
link_file(){
    local newfname="$2"
    local oldfname="$1"
    echo "$oldfname"
    if [ ! -f "$oldfname" ]
    then
        echo "ERROR:Cannot find $oldfname"
        exit 1
    fi
    if [ -h "$newfname" ]
    then
	if [ "$flag_lnfile" == "y" ]
        then
            if [ "$flag_lnfile_x" == "y" ]
            then
                rm $newfname
                ln -s $oldfname $newfname
            else
                return
            fi
        fi
        echo "$newfname already exists"
        echo "do you want to replace it ? (y/n)"
        read local x
	echo "do you want to force this option to all linked files (y/n)"
        read flag_lnfile
        if [ "$flag_lnfile" == "y" ]
        then
            flag_linfile_x=$x
        fi
        if [ "$x" == "y" ]
        then
            rm -rf $newfname
	    if [ -h "$newfname" ]
	    then
		echo "stupid"
		exit
	    fi
            ln -s $oldfname $newfname
        fi
    else
        ln -s $oldfname $newfname
    fi
}
create_dir "$subdir"
for subsubdir in "input" "output" "pdf" "tmp" "shower"
do
    create_dir "$subdir/$subsubdir"
    if [ "$subsubdir" == "input" ]
    then
	for cpfile in "process.inp" "default.inp" "user.inp" "seed.input" "decay_default.inp" "decay_user.inp" "decay_param_user.inp" "decay_param_default.inp" "shower_card_user.inp" "shower_card_default.inp" "fortran_compiler"
	do
	    copy_file "../$subsubdir/$cpfile" "$subdir/$subsubdir/$cpfile"
	done
	create_dir "$subdir/$subsubdir/lhapdf"
	for cpfile in "call_alphas" "call_strf_lhapdf"
	do
	    copy_file "../$subsubdir/lhapdf/$cpfile" "$subdir/$subsubdir/lhapdf/$cpfile"
	done
	create_dir "$subdir/$subsubdir/paths"
	for cpfile in $(ls -A "../$subsubdir/paths")
	do
	    copy_file "../$subsubdir/paths/$cpfile" "$subdir/$subsubdir/paths/$cpfile"
	done
    elif [ "$subsubdir" == "pdf" ]
    then
	for subsubsubdir in "cteq" "mrs" "gsdpdf"
	do
	    create_dir "$subdir/$subsubdir/$subsubsubdir"
	    if [ "$subsubsubdir" == "cteq" ]
	    then
		for cpfile in "cteq5l.tbl" "cteq5m.tbl" "cteq6d.tbl" "cteq6l.tbl" "cteq6l1.tbl" "cteq6m.tbl"
		do
		    copy_file "../$subsubdir/$subsubsubdir/$cpfile" "$subdir/$subsubdir/$subsubsubdir/$cpfile"
		done
	    elif [ "$subsubsubdir" == "mrs" ]
	    then
		for cpfile in "mrsb.dat" "mrse.dat" "mrst2002nlo.dat"
		do
		    copy_file "../$subsubdir/$subsubsubdir/$cpfile" "$subdir/$subsubdir/$subsubsubdir/$cpfile"
		done
	    elif [ "$subsubsubdir" == "gsdpdf" ]
	    then
		create_dir "$subdir/$subsubdir/$subsubsubdir/dpdfgrids"
		for cpfile in $(ls "../$subsubdir/$subsubsubdir/dpdfgrids")
		do
		    copy_file "../$subsubdir/$subsubsubdir/dpdfgrids/$cpfile" "$subdir/$subsubdir/$subsubsubdir/dpdfgrids/$cpfile"
		done
	    fi
	done
    elif [ "$subsubdir" == "shower" ]
    then
	for subsubsubdir in "HERWIG6" "HERWIGPP" "PYTHIA6" "PYTHIA8" "various"
	do
	    create_dir "$subdir/$subsubdir/$subsubsubdir"
	    for cpfile in $(ls "../$subsubdir/$subsubsubdir")
	    do
		copy_file "../$subsubdir/$subsubsubdir/$cpfile" "$subdir/$subsubdir/$subsubsubdir/$cpfile"
	    done
	done
	copy_file "../$subsubdir/RUN_HO_SHOWER" "$subdir/$subsubdir/RUN_HO_SHOWER"
    fi
done
if [ $# -ge 2 ]
then
    link_file "$hodir/bin/$addonfile" "$subdir/Helac-Onia"
else
    link_file "$hodir/bin/Helac-Onia" "$subdir/Helac-Onia" 
fi