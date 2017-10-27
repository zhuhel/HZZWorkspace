HULINDIR=/afs/cern.ch/user/x/xizhao/work/public/Results_2015-11-24_High/High
THEORYFILE=TheorySys.txt
POLYFILE=/afs/cern.ch/user/g/gcree/public/Analysis/PerEventFits/HZZWorkspace/polySys.txt
OUTPUTDIR=./output
MASSES="200 300 400 500 600 700 800 900 1000"
DUMMYMASSES=" 1200 1400 1600 1800 2000 2200 2400 2600 2800 3000"

rm -rf $OUTPUTDIR
mkdir $OUTPUTDIR


#
# renaming of txt to match out config files
#
for M in $MASSES;
do
 cp -f $HULINDIR/norm_ggH${M}_High.txt $OUTPUTDIR/ggH_${M}_norm.txt
 cp -f $HULINDIR/norm_VBFH${M}_High.txt $OUTPUTDIR/VBFH_${M}_norm.txt
 cp -f $HULINDIR/mean_ggH${M}_High.txt $OUTPUTDIR/ggH_${M}_shape_mean.txt
 cp -f $HULINDIR/mean_VBFH${M}_High.txt $OUTPUTDIR/VBFH_${M}_shape_mean.txt
done


#
# Hulins systematics only provided up to 1000GeV, so make dummy files above that (they don't enter result)
#
for M in $DUMMYMASSES;
do
  cp -f $OUTPUTDIR/ggH_1000_norm.txt $OUTPUTDIR/ggH_${M}_norm.txt
  cp -f $OUTPUTDIR/VBFH_1000_norm.txt $OUTPUTDIR/VBFH_${M}_norm.txt
  cp -f $OUTPUTDIR/ggH_1000_shape_mean.txt $OUTPUTDIR/ggH_${M}_shape_mean.txt
  cp -f $OUTPUTDIR/VBFH_1000_shape_mean.txt $OUTPUTDIR/VBFH_${M}_shape_mean.txt
done

MASSES=$MASSES$DUMMYMASSES


echo "doing theory systematics..."
#
# Now fill theory syst into the files
#
for M in $MASSES;
do
  LINE=`cat $THEORYFILE | grep "^$M "`

  SCALEDOWN=`echo $LINE | awk '{print $2}'`
  SCALEUP=`echo $LINE | awk '{print $3}'`
  PDFDOWN=`echo $LINE | awk '{print $4}'`
  PDFUP=`echo $LINE | awk '{print $5}'`

  TOINSERTQ="ATLAS_QCDscale_ggVV = $SCALEDOWN $SCALEUP"
  TOINSERTP="ATLAS_pdf_gg = $PDFDOWN $PDFUP"

  sed -i "s,],]\n${TOINSERTQ}\n${TOINSERTP}," $OUTPUTDIR/*_${M}_norm.txt

done 

#
# this business is adding polynomial acceptance systematic
#


#then fill in real deal
echo "now filling real acceptance sys.."
for M in $MASSES;
do

  for CAT in 4mu 2mu2e 4e;
  do
    for PROD in ggH VBFH;
    do
      POLDOWN=`cat $POLYFILE | grep $PROD | grep $CAT | grep " $M " | awk '{print $4}'`
      POLUP=`cat $POLYFILE | grep $PROD | grep $CAT | grep " $M " | awk '{print $5}'`
      if [ "$POLDOWN" == "" ]; then POLDOWN=1.0; fi
      if [ "$POLUP" == "" ]; then POLUP=1.0; fi

      echo "for $PROD $M $CAT want to fill $POLDOWN $POLUP" 
      
#sed -i -- "s,\(\[${PROD}_${CAT}_13TeV\]\),\1ATLAS_signal_acceptance_fit_${PROD}_${CAT} = $POLDOWN $POLUP," $OUTPUTDIR/${PROD}_${M}_norm.txt
      sed -i "s,\(\[ggF_${CAT}_13TeV\]\),\1\nATLAS_signal_acceptance_fit_${PROD}_${CAT} = $POLDOWN $POLUP," $OUTPUTDIR/${PROD}_${M}_norm.txt

    done
  done

done

