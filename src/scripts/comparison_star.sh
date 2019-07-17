#!/bin/bash

COMPDIR=$1
SRCDIR=$2
OUTDIR=$3

if [[ "$#" -ne 3 ]]; then
    echo "Usage: COMPDIR SRCDIR OUTDIR"
	exit 1
fi

mkdir -p ${OUTDIR}

for G in SRR3146803 SRR3146859 SRR3146914; do
	for S in s{1..5}; do
		for L in 100 500 1000 5000; do
			for T in circ_bsj circ_RF concordant discordant fusion keep many_hits no_hit OEA OEA2 orphan; do
				echo "bash ${COMPDIR}/ART/extract.sh $T $G $S $L ${SRCDIR}" > ${OUTDIR}/$G.$S.$L.$T.sh
				echo "bash ${COMPDIR}/STAR/extract.sh $T $G $S $L ${SRCDIR}" >> ${OUTDIR}/$G.$S.$L.$T.sh
				echo "bash ${COMPDIR}/COMPARISON/eval.sh $T $G $S $L ${SRCDIR}" >> ${OUTDIR}/$G.$S.$L.$T.sh
			done
		done
	done
done
