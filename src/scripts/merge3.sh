SRC=$1

echo "Merging results in working dir: $SRC"
mkdir -p ${SRC}/merged
DEST=${SRC}/merged

for M in concordant discordant circ_bsj circ_RF fusion keep OEA OEA2 orphan many_hits no_hit
	do
	echo "${M}"
	cat ${SRC}/*.${M}.pam > ${DEST}/${M}.pam
done

echo "Merged results:"
wc -l ${DEST}/*pam 
