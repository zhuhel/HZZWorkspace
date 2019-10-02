#!/bin/bash 
base=$HZZWSCODEDIR
files="fiducial_stat.root fiducial_sys.root"
totFiles="total_stat.root total_sys.root"
#put the folder where your resutls are 
subdirs="$base/xsec/results/latest"
isAsimov="no asimov"
isSyst="no sys" 


I=0
for f in $files ; do : 
    for a in $isAsimov ; do :  
	echo "#!/bin/bash" > job_"$I".sh 
	echo "cd $base"  >> job_"$I".sh
	echo "source setup.sh " >> job_"$I".sh
	if [[ $f == *sys* ]] ; then 
	    echo "fitFiducial $subdirs/$f $a sys" >> job_"$I".sh
	    chmod +x job_"$I".sh
	    bsub -q 1nd job_"$I".sh -o nll_"$I".out
	else
	    echo "fitFiducial $subdirs/$f $a no" >>job_"$I".sh
	    chmod +x job_"$I".sh
	    bsub -q 8nh job_"$I".sh -o nll_"$I".out
	    fi 
	I=$((I+1))
    done 
done 


for f in $totFiles ; do :
    for a in $isAsimov ; do :
        echo "#!/bin/bash" > job_"$I".sh
        echo "cd $base"  >> job_"$I".sh
        echo "source setup.sh " >> job_"$I".sh
        if [[ $f == *sys* ]] ; then 
            echo "fitFiducial $subdirs/$f $a sys total" >>job_"$I".sh
            chmod +x job_"$I".sh
            bsub -q 1nd job_"$I".sh -o nll_"$I".out
        else
            echo "fitFiducial $subdirs/$f $a no total" >>job_"$I".sh
            chmod +x job_"$I".sh
            bsub -q 8nh job_"$I".sh -o nll_"$I".out
        fi 
        I=$((I+1))
    done
done
