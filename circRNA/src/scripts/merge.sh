SRC=$1

echo "Merging results in working dir: $SRC"
mkdir -p ${SRC}/merged
DEST=${SRC}/merged

for N in 1 2 
	do
	echo "R${N}"
	cat ${SRC}/*ignore_R${N}.fastq > ${DEST}/ignore_R${N}.fastq
	for M in keep OEA orphan
		do
		cat ${SRC}/*_3.${M}_R${N}.fastq > ${DEST}/${M}_R${N}.fastq
	done
done

echo "Merged results:"
wc -l ${DEST}/*R1* | awk '{print $1/4"\t"$2}'
