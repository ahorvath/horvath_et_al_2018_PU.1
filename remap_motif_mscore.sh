#!/bin/bash

if [ "$#" -lt 4 ]; then
	echo "USAGE: $0 <PEAK> <SIZE> <GENOME> <MOTIF1> [<MOTIF2>] ..."
	exit -1
fi

peak=$1
size=$2
genome=$3
shift 3
motif=$@

motives=""; 
for motiffile in $motif; do 
	 motives=$motives"_`basename $motiffile .motif`"; 
done

motif_count=$#
last_col=`expr 22 + ${motif_count}`

sample=`basename $peak .bed`${motives}_${genome}_s_${size}
mkdir -p ${sample}

annotatePeaks.pl $peak ${genome} -size $size -m $motif -mscore -mbed ${sample}/${sample}_mbed.bed | tee ${sample}/${sample}_w_motives.txt | sed '1d' | cut -f2,3,4,22-${last_col} | sortBed | cat <(echo -e "Chr\tStart\tEnd$motives" |  sed 's/_/\t/g') - > ${sample}/${sample}_w_motives.tsv 2> ${sample}/${sample}.log
sed '1d' ${sample}/${sample}_w_motives.tsv > ${sample}/${sample}_w_motives.bed
totalpeak=`wc -l $peak | cut -d" " -f1`;
motivespeak=`wc -l ${sample}/${sample}_w_motives.bed | cut -d" " -f1`;
rate=`echo "scale=4;$motivespeak/$totalpeak" | bc`;
echo -e "$totalpeak\t$motivespeak\t$rate";
