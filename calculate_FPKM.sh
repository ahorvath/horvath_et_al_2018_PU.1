if [[ $# -le 1 ]]; then
       echo -e "USAGE: $0 <BAM> <BED> \n"
       echo -e "\tEXAMPLE: $0 sample.bam sample.bed > sample_FPKM.tsv"
       exit 1
fi

SAMPLE=$1
BED=$2

mkdir -p {stats,counts};
>&2 echo "SAMPLE $SAMPLE"; 
sample=`basename $SAMPLE .bam`;
bamtools stats -in $SAMPLE > stats/${sample}.bamstats;

coverageBed -counts -abam $SAMPLE -b $BED > counts/${sample}_counts.bed;
        
>&2 echo "Sample name $sample";
reads=`grep Mapped stats/${sample}.bamstats | awk '{print $3}'`
>&2 echo "Read count: $reads"
awk -v reads=$reads -F"\t" '{OFS="\t"; print $1,$2,$3,($NF*1000000000)/((reads+1)*($3-$2))}' counts/${sample}_counts.bed | sortBed 

