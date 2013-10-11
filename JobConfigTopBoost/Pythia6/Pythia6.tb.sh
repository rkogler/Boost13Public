#!/bin/bash -l

TMPDIR=${_CONDOR_SCRATCH_DIR:-/tmp}
JOB_ID=$1
TUNE=$2
NEV=10000

RIVET_PREFIX=${HOME}/Rivet
RIVET_ANALYSIS_DIR=${RIVET_PREFIX}/Analysis/rivet-cflow
RIVET_ANALYSIS="MC_CFLOW_TTBAR -a MC_SUBSTRUCTURE"
OUTDIR=${RIVET_ANALYSIS_DIR}/HighStatsAida/Pythia6
mkdir -p ${OUTDIR}

export LHAPATH=${RIVET_PREFIX}/local/share/lhapdf/PDFsets
export AGILE_GEN_PATH=${RIVET_PREFIX}/local/generators 
source ${RIVET_PREFIX}/rivetenv.sh
source ${RIVET_PREFIX}/agileenv.sh
STATCODE=0

#remake files with bad binning in jcharge. all files exist.
if [  ! -s ${OUTDIR}/Pythia6.tb._ptune${TUNE}_part${JOB_ID}.aida -o \
      ${OUTDIR}/Pythia6.tb._ptune${TUNE}_part${JOB_ID}.aida -ot ${OUTDIR}/Pythia6.tb._ntune354_part98.aida ]; then 
    SEED_PARAM=`echo "MRPY(1)=1270${JOB_ID}"`
    cd ${RIVET_ANALYSIS_DIR}
    FIFONAME=${TMPDIR}/pythia6pb_fifo_${TUNE}_${JOB_ID}.hepmc
    [ -e ${FIFONAME} ] && rm ${FIFONAME} #cleanup old jobs
    mkfifo ${FIFONAME}
    agile-runmc Pythia6:426 --beams=LHC:7000 -n $NEV \
	-P fpythia-ttbar-poslepton.params \
	-p PYTUNE=${TUNE}  -p "CKIN(3)=200" -p ${SEED_PARAM} -o ${FIFONAME} &
    sleep 1 #let things settle after starting AGILe
    rivet --analysis-path=${RIVET_ANALYSIS_DIR} -a ${RIVET_ANALYSIS} \
	--histo-file=${OUTDIR}/Pythia6.tb._ptune${TUNE}_part${JOB_ID}.aida ${FIFONAME}
    STATCODE=$? 
    [ -e ${FIFONAME} ] && rm ${FIFONAME} #cleanup this job
fi

if [  ! -s ${OUTDIR}/Pythia6.tb._ntune${TUNE}_part${JOB_ID}.aida -o \
      ${OUTDIR}/Pythia6.tb._ntune${TUNE}_part${JOB_ID}.aida -ot ${OUTDIR}/Pythia6.tb._ntune354_part98.aida  ]; then 
    SEED_PARAM=`echo "MRPY(1)=1240${JOB_ID}"`
    FIFONAME=${TMPDIR}/pythia6nb_fifo_${TUNE}_${JOB_ID}.hepmc
    [ -e ${FIFONAME} ] && rm ${FIFONAME} #cleanup old jobs
    mkfifo ${FIFONAME}
    agile-runmc Pythia6:426 --beams=LHC:7000 -n $NEV \
	-P fpythia-ttbar-neglepton.params \
	-p PYTUNE=${TUNE} -p "CKIN(3)=200" -p ${SEED_PARAM} -o ${FIFONAME} &
    sleep 1 #let things settle after starting AGILe
    rivet --analysis-path=${RIVET_ANALYSIS_DIR} -a ${RIVET_ANALYSIS} \
	--histo-file=${OUTDIR}/Pythia6.tb._ntune${TUNE}_part${JOB_ID}.aida ${FIFONAME}
    STATCODE=$? 
    [ -e ${FIFONAME} ] && rm ${FIFONAME} #cleanup this job
fi

exit $STATCODE
return 0 #does this prevent the end-of-job-hang?

