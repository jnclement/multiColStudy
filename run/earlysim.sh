#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
STARTN=$(( $2 * 4 ))
ISSIM=$4
RN=$5
for i in {0..3}; do
    SUBDIR=$(( $STARTN + $i ))
    UPLN=$(( $SUBDIR + 1 ))
    mkdir -p $SUBDIR
    mkdir -p $SUBDIR\_multicol
    mkdir -p lists
    mkdir -p ./output/smg
    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/${SUBDIR}/
    mkdir -p /sphenix/tg/tg01/jets/jocl/multiCol/${RN}/
    mkdir -p ./dsts/$SUBDIR
    ls -larth
    cp -r /sphenix/user/jocl/projects/multiColStudy/run/run_earlydata.C .
    if [ $ISSIM -ne 0 ]; then
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/g4hits.list ./lists/g4hits.list
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_truth_jet.list ./lists/dst_truth_jet.list
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_calo_cluster.list ./lists/dst_calo_cluster.list
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_global.list ./lists/dst_global.list
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_mbd_epd.list ./lists/dst_mbd_epd.list
    else
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_calofitting_run2pp-000${RN}.list ./lists/dst_calo_cluster.list
	cp -r /sphenix/user/jocl/projects/multiColStudy/run/lists/dst_triggered_event_run2pp-000${RN}.list ./lists/dst_global.list
    fi
    G4HITSF=`sed -n "${UPLN}"p ./lists/g4hits.list`
    CALOCLF=`sed -n "${UPLN}"p ./lists/dst_calo_cluster.list`
    SEGEND=`echo $CALOCLF | awk -F"-" '{print $3}'`
    GLOBALF=""
    if [ $ISSIM -eq 0 ]; then
	GLOBALF="DST_TRIGGERED_EVENT_run2pp_ana462_nocdbtag_v001-000${RN}-${SEGEND}"
    else
	GLOBALF=`sed -n "${UPLN}"p ./lists/dst_global.list`
    fi
    TRTHJET=`sed -n "${UPLN}"p ./lists/dst_truth_jet.list`
    DMBDEPD=`sed -n "${UPLN}"p ./lists/dst_mbd_epd.list`
    FULLTRTH=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${TRTHJET}';"`
    FULLMBEP=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${DMBDEPD}';"`
    FULLCALO=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${CALOCLF}';"`
    FULLG4HT=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${G4HITSF}';"`
    FULLGLOB=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${GLOBALF}';"`
    echo $CALOCLF
    #echo $GLOBALF
    #getinputfiles.pl $GLOBALF
    #getinputfiles.pl $CALOCLF
    #getinputfiles.pl $TRTHJET
    #getinputfiles.pl $DMBDEPD
    #getinputfiles.pl $G4HITSF
    cp $FULLCALO .
    
    if [ $ISSIM -ne 0 ]; then
	cp $FULLGLOB .
	cp $FULLTRTH .
	cp $FULLMBEP .
	cp $FULLG4HT .
    fi
    #cp -r $G4HITSF ./dsts/$2/g4hits_${2}.root
    echo ""
    echo "" 
    ls
    echo ""
    echo ""
    mv $CALOCLF ./dsts/$SUBDIR/calo_cluster_${SUBDIR}.root

    if [ $ISSIM -ne 0 ]; then
	mv $GLOBALF ./dsts/$SUBDIR/global_${SUBDIR}.root
	mv $G4HITSF ./dsts/$SUBDIR/g4hits_${SUBDIR}.root
	mv $TRTHJET ./dsts/$SUBDIR/truth_jet_${SUBDIR}.root
	mv $DMBDEPD ./dsts/$SUBDIR/mbd_epd_${SUBDIR}.root
    fi
    ls ./dsts/$SUBDIR
    #cp -r $TRTHJET ./dsts/$SUBDIR/truth_jet_${SUBDIR}.root
    root -b -q 'run_earlydata.C("'${1}'",'${SUBDIR}',0,'${3}',".",'${ISSIM}','${RN}')'
    ls
    echo " "
    ls $SUBDIR
    echo " "
    ls $SUBDIR\_multicol -larth
    echo "copying"
    if [ $ISSIM -ne 0 ]; then
	cp -r ./${SUBDIR}/* /sphenix/tg/tg01/jets/jocl/multiCol/${SUBDIR}/
	cp -r ./${SUBDIR}_multicol/* /sphenix/tg/tg01/jets/jocl/multiCol/${SUBDIR}/
    else
	cp -r ./${SUBDIR}/* /sphenix/tg/tg01/jets/jocl/multiCol/${RN}/
	cp -r ./${SUBDIR}_multicol/* /sphenix/tg/tg01/jets/jocl/multiCol/${RN}/
    fi
    echo "done with loop ${SUBDIR}"
done

