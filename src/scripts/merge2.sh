SRC=$1

echo "Merging results in working dir: $SRC"
mkdir -p ${SRC}/merged
DEST=${SRC}/merged

CONCRD=concordant_cdna

for N in 1 2 
	do
	echo "R${N}"
	cat ${SRC}/*${CONCRD}_R${N}.fastq > ${DEST}/${CONCRD}_R${N}.fastq
	for M in concordant_gen discordant circ_bsj circ_RF fusion keep OEA OEA2 orphan many_hits no_hit
		do
		echo "${M}"
		cat ${SRC}/*_3.${M}_R${N}.fastq > ${DEST}/${M}_R${N}.fastq
	done
done

cat  ${SRC}/*.junc > ${DEST}/all.junc

echo "Merged results:"
wc -l ${DEST}/*R1* | awk '{print $1/4"\t"$2}'
