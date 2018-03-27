#!/bin/sh

# Author: Endre Barta & Gergely Nagy
# v0.1 2010/12/25
# v1.0 2011/03/08
# v1.1 2011/05/09 support for EBI ENA SRA fastq.gz ftp location added
# v1.2 2011/05/21 alternative bwa command for trimming barcodes added
# v1.3 2011/07/18 changing some motif finding parameters for homer2 (can use more processor)
# v1.4 2011/10/21 introducing the third field in the list file 
#	HOMER2's makeTagDirectory can accept now BAM format alignment files
# v1.4.1 2011/12/04 generating larger and strand specific bedgraph file
#	bed and bedgraph files go to the $basedir/bed  and bedgraph directories
# v1.5 2012/02/06 processing GRO-seq results as well using HOMER
#	will generate transcripts, double and single peaks and transcripts counts
#	I didn't update the ChIPseeqer code, it is now commented out
# v1.6 2012/06/04 small changes, bugfixes like writing out individual tmp files and creating bigbasedir if needed
# v1.7 2012/07/13 upgrading to MACS2 and some other small changes
# v1.8 2012/10/01 HOMER denovo is now using the top 1000 peaks centered 
# v1.9mm10 default mouse refrence genome is now mm10

# getting a list of ChIP-seq experiments from the list file
# first field is the short name of the experiment (e.g. mm_L1_PPARg1_d7)
# The first two characters must be hs or mm. If not, hs is assumed
# second field is the ftp path to the SRA-LITE (SRX) directory where the sralite (*.sra) format raw data are
# the second field can also be the location for the EBI ENA SRR directory like:
#	ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR039/SRR039818
# it can be omitted or just put a tab instead if the fastq or BAM files already are in place
# the third field can have the following values:
#	factor	# it is the default (if nothing is written, or anything else that listed
#	histone	# it is a histone modification
#	control # it indicates that this chip has to be used for control in the following experiments
#	groseq # it is a GRO-seq experiment MACS and ChIPseeqer processing will be omitted

# the second argument is the base directory for the analysis

# the third argument is the base directory for the large files (.sra, .fastq.gz, .sai) needed for the analysis
# this can be the same as the base directory

# parameters that can be changed:
procno=3 # number of processors (cores) to be used by BWA for aligning and meme_p for de novo motif finding
	 # number of available processors (cores) minus 1 is a wise choice (if other programs like HOMER motif finding are not running in the background)
#motif1=rxr # for a list of available motifs see the data/knownTFs/motifs in the homer root directory
#motif2=pu1 # motifs to annotate with HOMER
#motif2=stat6 # motifs to annotate with HOMER
#cmotif1=PPARG-RXRA #for a list of available motifs see $CHIPSEEQERDIR/DATA/JASPAR_MATRICES
#cmotif2=Cebpa # motifs to annotate with ChIPseeqer
motiflength="10,13,16" #length (or lengths of motifs) to look for (HOMER) can be in the following form: 8,10,13
#minmotifl=10 # minimum length of motifs to find for MEME
#maxmotifl=13 # maximum length
#nmotif=3 #number of motifs to find for MEME
#--- homer2 does not use this ---# depth=med # for HOMER motif finding [low|med|high|allnight] (time spent on local optimization default: med)
depth=`echo $motiflength | sed "s/,/_/g"`
indexdir=/molbio/flatfiles/genomes/index # directory for BWA indexed genome files
#phastdir=/molbio/flatfiles/genomes/hs/hg18/phastCons44way/placentalMammals # needed only for ChIPseeqer human conservation analysis
barcode_length=8 # I didn't check whether the 0 should work (if yes, that should be the default)
		# Otherwise, trimming the first (in most cases the poorest quality) base does not hurt
style=factor	# for homer style can be factor or histone
iscontrol="" # by default there is no control (IG) experiment
macscontrol=""
control=""
#fragLength=" -fragLength auto "
fragLength=" -fragLength 150 " # for HOMER's makeUCSCfile it is better not to let it to guess but keep it standard
listfile=$1 # first argument is the list file
basedir=$2 # second argument is the base directory for the analysis
if [ $3 ] ; then
	bigbasedir=$3 # third argument for storing the large analysis files
else
	bigbasedir=$basedir
fi

if [ ! -d $bigbasedir ] ; then 
	/bin/mkdir -p $bigbasedir
fi

logbasedir=${basedir}/logs
if [ ! -d $logbasedir ] ; then
	mkdir $logbasedir
fi

if [ ! -d ${basedir}/bed ] ; then mkdir ${basedir}/bed ; fi

for i in SRA fastq sai ; do # making directories for the large files if not existed
	if [ ! -e ${bigbasedir}/$i ] ; then
	       mkdir ${bigbasedir}/$i
       fi
done

for i in `grep "^[^#]" $listfile | sed 's/ /;/g' | sed "s/\t/;/g"` ; do
	name=`echo $i | awk -F";" '{print $1}'`
	sp=`echo $name | cut -c 1-2`
	homerstyle=`echo $i | awk -F";" '{print $3}'`
	if [ "$homerstyle" = "control" ] ; then #test the third field in the listfile and set the style parameter for the HOMER's findPeaks
		#style=histone	#let it be factor at the beginning and later always the value of the previous sample
		iscontrol=1
	elif [ "$homerstyle" = "histone" ] ; then
		style=histone
	else
		style=factor
	fi
	if [ "$sp" = "mm" ] ; then
		species="mm"
		genome="mm10"
		spec="mouse"
		uchr="chrX chrY"
	elif [ "$sp" = "hs" ] ; then
		species="hs"
		genome="hg19"
		spec="human"
		uchr="chr20 chr21 chr22 chrX chrY"
	else
                species="mm"		#The default is mouse (when the form is not hs_ or mm_)
                genome="mm10"
                spec="mouse"
                uchr="chrX chrY"
	fi

	#logdir=$logbasedir
	logdir=${logbasedir}/$name # Using this more structured option
	if [ ! -d $logdir ] ; then
		mkdir $logdir
	fi
	ftpsite=`echo $i | awk -F";" '{print $2}'`
	if [ -d ${basedir}/${name} ] ; then
		#/bin/mv ${basedir}/${name} ${basedir}/${name}.old
		echo "$name already processed at some extent"
		cd ${basedir}/${name}
		if [ ! -L logs ] ; then
			ln -s ${logdir} logs
		fi
	else
		mkdir ${basedir}/${name}
		cd ${basedir}/${name}
		ln -s ${logdir} logs
	fi
	if [ ! -d ${bigbasedir}/SRA/$name ] ; then
		mkdir ${bigbasedir}/SRA/$name
	fi
	if [ ! -L ${basedir}/${name}/SRA ] ; then
		ln -s ${bigbasedir}/SRA/$name SRA
	fi
	for d in fastq sai ; do #we are at ${basedir}/${name}
		if [ ! -e ${basedir}/${name}/$d ] ; then #basedir not equal bigbasedir and first analysis of $name
			ln -s ${bigbasedir}/$d $d
		fi
	done
	#for d in bam macs homer meme bed chipseeqer ; do #no need for chipseeqer, meme and bed
	for d in bam macs homer ; do #no need for chipseeqer
		if [ ! -d ${basedir}/${name}/$d ] ; then
			mkdir $d 
		fi
	done
	if [ ! -d ${basedir}/bedgraph ] ; then
		mkdir ${basedir}/bedgraph
	fi
	if [ ! -d ${basedir}/fasta ] ; then
		mkdir ${basedir}/fasta # for the peak and beackground sequence files
	fi
	if [ -d ${basedir}/${name}/fastq/$name ] ; then
		/bin/rm -rf ${basedir}/${name}/fastq/$name # this directory holds only the fastq files generated by sra-lite
	fi
	mkdir ${basedir}/${name}/fastq/$name
	# finished making directory structure

	cd ${basedir}/${name}/SRA
	echo "listfile: $i"
	echo "name=$name sp=$sp style=$style icontrol=$iscontrol genome=$genome ftpsite=$ftpsite"
	echo "starting processing $name at `date` ..." > ${logdir}/${name}_ChIP-seq_analyze.log
	if [ ! -f ${basedir}/${name}/fastq/${name}.fastq.gz -a ! -f ${basedir}/${name}/bam/${name}.bam ] ; then
		if [ ! -n "$ftpsite" ] ; then #no ftpsite, no fastq file, no BAM file = ERROR!
			echo "Something wrong, at least one of the ftpsite or the fastq file or the BAM file must be present" >> ${logdir}/${name}_ChIP-seq_analyze.log
			echo "name=$name sp=$sp style=$style icontrol=$iscontrol genome=$genome ftpsite=$ftpsite"
			exit
		fi
		wget -o ${logdir}/${name}-wget.log --dns-timeout=100 -r -m -nH -nd ${ftpsite}/*
		echo "finished wget $ftpsite" >> ${logdir}/${name}_ChIP-seq_analyze.log
		echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
		# finished wget, starting fastq-dump

		if ls -l | grep -q ".fastq.gz$" ; then
			# if the transferred files are from EBI, move them to the fastq directory:
			for i in *.fastq.gz ; do gunzip -c $i ; done |
				gzip -c > ${basedir}/${name}/fastq/${name}.fastq.gz # this needs if there is more than one run
				/bin/rm *.fastq.gz
			# done
		else	
			echo "Starting converting SRA files to fastq format" >> ${logdir}/${name}_ChIP-seq_analyze.log
			cd ${basedir}/${name}/fastq/${name}
			for sra in ${basedir}/${name}/SRA/*.sra ; do
				echo "fastq-dumping $sra for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log 
				fastq-dump $sra > ${logdir}/${name}-fastq-dump.log 2> ${logdir}/${name}-fastq-dump.err
			done
		fi
		/bin/rm -rf ${bigbasedir}/SRA/${name}  # remove SRA directory compeletly
		/bin/rm -f ${basedir}/${name}/SRA  # remove SRA directory compeletly
		#/bin/rm 
		echo "finished fastq-dumping for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
		echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
		# finished fastq-dump starting indexing (generating sai)
		echo "Starting gzipping fastq files" >> ${logdir}/${name}_ChIP-seq_analyze.log
		cd ${basedir}/${name}/fastq
		cat ${basedir}/${name}/fastq/${name}/*.fastq | gzip -c > ${basedir}/${name}/fastq/${name}.fastq.gz
		/bin/rm -rf ${basedir}/${name}/fastq/$name
	fi
	/bin/rm -f ${basedir}/${name}/SRA #ususally wo no not need this
	if [ ! -f  ${basedir}/${name}/bam/${name}.bam ] ; then 
		echo "Indexing to sai for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
		cd ${basedir}/${name}/sai
		#( nice -n 20 bwa aln -t $procno  ${indexdir}/$genome ${basedir}/${name}/fastq/${name}.fastq.gz > ${name}.sai ) >  ${logdir}/${name}-sai.log 2> ${logdir}/${name}-sai.err
		( nice -n 20 bwa aln -t $procno -B $barcode_length ${indexdir}/$genome ${basedir}/${name}/fastq/${name}.fastq.gz > ${name}.sai ) >  ${logdir}/${name}-sai.log 2> ${logdir}/${name}-sai.err
		echo "finished aligning fastq for $genome for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
		echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
		# finished fastq alignment starting generating sam files
		echo "Starting converting sai files into sam/bam format" >> ${logdir}/${name}_ChIP-seq_analyze.log
		cd ${basedir}/${name}/bam
		bwa samse ${indexdir}/$genome ${basedir}/${name}/sai/${name}.sai ${basedir}/${name}/fastq/${name}.fastq.gz |
		samtools view -buS -t ${indexdir}/$genome.fai - |
		samtools sort -m 1000000000 - $name
		samtools index ${name}.bam ${name}.bai
		/bin/rm -f ${basedir}/${name}/sai/${name}.sai
		echo "finished converting sai to sam and bam finally on pipeline for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
		echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
	fi
	/bin/rm -f ${basedir}/${name}/sai # remove the fastq directory
	# finished alignment of the reads to the genome, starting MACS processing
	echo "Starting MACS peak finding process" >> ${logdir}/${name}_ChIP-seq_analyze.log
	cd ${basedir}/${name}/macs
	if [ -f ${name}_macs_peaks.xls ] ; then
		echo "MACS peak calling for $name already done" >> ${logdir}/${name}_ChIP-seq_analyze.log
	else
		#( macs2 callpeak -f BAM -n ${name}_macs -g $species -B --SPMR -t ${basedir}/${name}/bam/${name}.bam $macscontrol > ${logdir}/${name}-macs.log 2> ${logdir}/${name}-macs.err ; /bin/cp *.bed ${basedir}/bed/ ) &
		#macs2 callpeak -f BAM -n ${name}_macs -g $species -B --SPMR -t ${basedir}/${name}/bam/${name}.bam $macscontrol > ${logdir}/${name}-macs.log 2> ${logdir}/${name}-macs.err
		macs2 callpeak -f BAM -n ${name}_macs2 -g $species -t ${basedir}/${name}/bam/${name}.bam $macscontrol > ${logdir}/${name}-macs.log 2> ${logdir}/${name}-macs.err
		intersectBed -v -a ${name}_macs2_peaks.bed -b /molbio/flatfiles/${genome}_encode_blacklist.bed > ${basedir}/bed/${name}_macs_peaks.bed
		intersectBed -v -a ${name}_macs2_summits.bed -b /molbio/flatfiles/${genome}_encode_blacklist.bed > ${basedir}/bed/${name}_macs_summits.bed
		#sort -k5,5gr ${basedir}/bed/${name}_macs_summits.bed | head -1000 > ${basedir}/bed/${name}_macs_top.bed
		#/bin/cp *.bed ${basedir}/bed/
		#TODO correct it to ${name}_macs
		#/bin/cp ${name}_macs_peaks.encodePeak ${basedir}/bed/${name}_macs_narrow.bed
		#/bin/mv ${name}_macs_control_lambda.bdg ${name}_macs_control_lambda.bedgraph
		#gzip ${name}_macs_control_lambda.bedgraph
		#/bin/mv ${name}_macs_control_lambda.bedgraph.gz ${basedir}/bedgraph/
		#/bin/mv ${name}_macs_treat_pileup.bdg ${name}_macs_treat_pileup.bedgraph
		#gzip ${name}_macs_treat_pileup.bedgraph
		#/bin/mv ${name}_macs_treat_pileup.bedgraph.gz ${basedir}/bedgraph/
		echo "finished MACS peak calling for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
	fi
        echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
	# finished MACS processing, starting HOMER
	echo "Starting HOMER processing" >> ${logdir}/${name}_ChIP-seq_analyze.log
	cd ${basedir}/${name}/homer
	# first, need to convert the BAM file to bed format #--- obsolete, if samtools installed the program can use bam alignment files
	# putting the bed file to the fastq directory
	if [ ! -d ${name} ] ; then
	# if [ ! -f ${basedir}/${name}/fastq/${name}.bed ] ; then 
		# No need it any more # bamToBed -i ${basedir}/${name}/bam/${name}.bam -tag | awk -F"\t" '{ print $4,$1,$2,$3,$6}' | uniq -f1 -c | awk '{OFS="\t" ; print $3,$4,$5,$2,$1,$6}' > ${basedir}/${name}/fastq/${name}.bed
		makeTagDirectory ${name} ${basedir}/${name}/bam/${name}.bam -genome $genome -name ${name} -format sam default -checkGC > ${logdir}/${name}-maketagdirectory.log 2> ${logdir}/${name}-maketagdirectory.err
		#checkTagBias.pl $name $genome > ${logdir}/${name}-checkTagBias.log 2> ${logdir}/${name}-checkTagBias.err # = -checkGC above
		tagnum=`awk '/^genome=/ {print $3}' ${name}/tagInfo.txt | awk -F"." '{print $1}'`
		if [ "$homerstyle" = "groseq" ] ; then
			fragLength=" -fragLength 120 " #adjust it into the GRO-seq fragment length
		fi
		makeUCSCfile ${name} $control -name ${name} -o ${basedir}/bedgraph/${name}_big.bedgraph -fsize 1e50 -norm $fragLength > ${logdir}/${name}-makeUCSCfile_orig.log 2> ${logdir}/${name}-makeUCSCfile_orig.err
		igvtools toTDF ${basedir}/bedgraph/${name}_big.bedgraph.gz ${basedir}/bedgraph/${name}_big.bedgraph.tdf ${genome} > ${logdir}/${name}-toTDF.log
		#makeUCSCfile ${name} $control -name ${name}_small -o ${basedir}/bedgraph/${name}_small.bedgraph -norm $fragLength > ${logdir}/${name}-makeUCSCfile_orig_small.log 2> ${logdir}/${name}-makeUCSCfile_orig_small.err
                makeUCSCfile ${name} $control -fsize 1e50 -strand '+' $fragLength -norm -name ${name}_norm -style chipseq 2> ${logdir}/${name}-makeUCSCfile_forward.err | grep -v "^chrM" > /tmp/${name}bedgraph
		makeUCSCfile ${name} $control -fsize 1e50 -strand '-' -neg $fragLength -norm -name $name -style chipseq 2> ${logdir}/${name}-makeUCSCfile_reverse.err | grep -v "^track" | grep -v "^chrM" >> /tmp/${name}bedgraph
                head -1 /tmp/${name}bedgraph | sed "s/+ strand//" > ${basedir}/bedgraph/${name}.bedgraph
                for j in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 $uchr ; do #TODO make it works for human (chrs 20-22)
			grep "^${j}[^0-9]" /tmp/${name}bedgraph | sort -n -k2  >> ${basedir}/bedgraph/${name}.bedgraph
		done
		gzip ${basedir}/bedgraph/${name}.bedgraph
		if [ "$homerstyle" = "groseq" ] ; then
			# for GRO-seq we would need an un-normalized bedgraph as well
			makeUCSCfile ${name} -fsize 1e50 -strand '+' $fragLength -name ${name} -noadj -style chipseq 2> ${logdir}/${name}-makeUCSCfile_forward.err | grep -v "^chrM" > /tmp/${name}bedgraph
			makeUCSCfile ${name} -fsize 1e50 -strand '-' -neg $fragLength -name $name -noadj -style chipseq 2> ${logdir}/${name}-makeUCSCfile_reverse.err | grep -v "^track" | grep -v "^chrM" >> /tmp/${name}bedgraph
			head -1 /tmp/${name}bedgraph | sed "s/+ strand//" > ${basedir}/bedgraph/${name}_nonorm.bedgraph
			for j in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 $uchr ; do
				grep "^${j}[^0-9]" /tmp/${name}bedgraph | sort -n -k2  >> ${basedir}/bedgraph/${name}_nonorm.bedgraph
			done
			gzip ${basedir}/bedgraph/${name}_nonorm.bedgraph
			# for hg19 and mm9 (no uniqmap for mm10) findPeaks $name -style groseq -tssSize 200 -minBodySize 300 -maxBodySize 3000 -tssFold 5 -bodyFold 2 -endFold 10 -fragLength 35 -uniqmap /usr/local/molbio/homer/uniqmap/${genome}-uniqmap -o auto > ${logdir}/${name}-findPeaks.log 2> ${logdir}/${name}-findPeaks.err
			findPeaks $name -style groseq -tssSize 200 -minBodySize 300 -maxBodySize 3000 -tssFold 5 -bodyFold 2 -endFold 10 -fragLength 35 -o auto > ${logdir}/${name}-findPeaks.log 2> ${logdir}/${name}-findPeaks.err
			pos2bed.pl ${name}/transcripts.txt > ${basedir}/bed/${name}-homertranscripts.bed
			analyzeRNA.pl rna ${genome} -d $name -rpkm -count genes -strand + -normMatrix 1e6 > ${name}_rpkm_gene.txt 2> ${logdir}/${name}_rpkm_gene.err
			analyzeRNA.pl rna ${genome} -d $name -rpkm -count genes -normMatrix 1e6 -pausing 300 > ${name}_rpkm_pausing.txt 2> ${logdir}/${name}_rpkm_pausing.err
		elif [ "$homerstyle" = "histone" ] ; then
			# findPeaks parameter according to Chris's recommendations
			findPeaks $name $control -size 1000 -minDist 2500 -o ${name}_homerpeaks.txt >  ${logdir}/${name}-findPeaks.log 2> ${logdir}/${name}-findPeaks.err #for denovo motif finding
			findPeaks $name $control -region -size 1000 -minDist 2500 -o ${name}_region_homerpeaks.txt >  ${logdir}/${name}_region-findPeaks.log 2> ${logdir}/${name}_region-findPeaks.err
			pos2bed.pl ${name}_region_homerpeaks.txt > ${basedir}/bed/${name}_region_homerpeaks.bed
			#annotatePeaks.pl ${name}_region_homerpeaks.txt $genome -size 200 -go ${name}_region_$depth -genomeOntology ${name}_region_homermotifs_$depth \
                        #        -matrix ${name}_region_homer_motifsmatrix.txt -mbed ${name}_region_homer_motifs.bed \
                        #        -m /usr/local/molbio/homer/data/knownTFs/motifs/${motif1}.motif /usr/local/molbio/homer/data/knownTFs/motifs/${motif2}.motif \
                        #        > ${name}_region-homermotifsannot.txt 2> ${logdir}/${name}_region-annotatePeaks.err
		else
			findPeaks $name -center $control -style $style -o ${name}_homerpeaks.txt >  ${logdir}/${name}-findPeaks.log 2> ${logdir}/${name}-findPeaks.err
		fi
		if [ -f ${name}_homerpeaks.txt ] ; then
			findPeaks.pl $name $control -auto -ucsc ${name}_both -genome $genome -${spec} -dual -name $name -strand > ${logdir}/${name}-findPeakspl.log 2>${logdir}/${name}-findPeakspl.err
			#annotatePeaks.pl ${name}_homerpeaks.txt $genome -size 200 -go  ${name}_homermotifs_$depth -genomeOntology ${name}_homermotifs_$depth \
			#	-matrix ${name}_homer_motifsmatrix.txt -mbed ${name}_homer_motifs.bed \
			#	-m /usr/local/molbio/homer/data/knownTFs/motifs/${motif1}.motif /usr/local/molbio/homer/data/knownTFs/motifs/${motif2}.motif \
			#	> ${name}_homermotifsannot.txt 2> ${logdir}/${name}-annotatePeaks.err
                        #annotatePeaks.pl ${basedir}/${name}/macs/${name}_macs_peaks.bed $genome -size 200 -go ${name}_macs_homermotifs_$depth -genomeOntology ${name}_macs_homermotifs_$depth \
                        #        -matrix ${name}_macs_homer_motifsmatrix.txt -mbed ${name}_macs_homer_motifs.bed \
                        #        -m /usr/local/molbio/homer/data/knownTFs/motifs/${motif1}.motif /usr/local/molbio/homer/data/knownTFs/motifs/${motif2}.motif \
                        #        > ${name}_macs-homermotifsannot.txt 2> ${logdir}/${name}_macs-annotatePeaks.err
			pos2bed.pl ${name}_homerpeaks.txt | intersectBed -v -a stdin -b /molbio/flatfiles/${genome}_encode_blacklist.bed > ${basedir}/bed/${name}-homerpeaks.bed
			sort -k5,5gr ${basedir}/bed/${name}-homerpeaks.bed | head -1000 > ${basedir}/bed/${name}-homer_top.bed
		fi
		#if [ "$homerstyle" = "histone" ] ; then
			#nice -n 20 findMotifsGenome.pl ${name}_homerpeaks.txt ${genome} ${name}_homermotifs_$depth -size given -mask -len $motiflength -p $procno -dumpFasta -bits -preparse -homer2 > ${logdir}/${name}-findMotifsGenome_${depth}.log 2> ${logdir}/${name}-findMotifsGenome_${depth}.err
		#el
		if [ "$homerstyle" = "factor" -o "$homerstyle" = "control" ] ; then
			nice -n 20 findMotifsGenome.pl ${basedir}/bed/${name}-homer_top.bed ${genome} ${name}_homermotifs_$depth -mask -size 100 -len $motiflength -p $procno -dumpFasta -bits -preparse -homer2 > ${logdir}/${name}-findMotifsGenome_${depth}.log 2> ${logdir}/${name}-findMotifsGenome_${depth}.err
			/bin/mv ${name}_homermotifs_${depth}/target.fa ${basedir}/fasta/${name}_homerpeaks.dfasta
			gzip ${basedir}/fasta/${name}_homerpeaks.dfasta
			/bin/mv ${name}_homermotifs_${depth}/background.fa ${basedir}/fasta/${name}_homerbackground.dfasta
			gzip ${basedir}/fasta/${name}_homerbackground.dfasta
			#nice -n 20 findMotifsGenome.pl ${basedir}/${name}/macs/${name}_macs_sum100top.bed ${genome} ${name}_macs_homermotifs_$depth -mask -size given -len $motiflength -p $procno -dumpFasta -bits -preparse -homer2 > ${logdir}/${name}_macs-findMotifsGenome_${depth}.log 2> ${logdir}/${name}_macs-findMotifsGenome_${depth}.err # doing the denovo motif finding with the MACS peaks as well 
                	#/bin/mv ${name}_macs_homermotifs_${depth}/target.fa ${basedir}/fasta/${name}_macspeaks.dfasta
                	#gzip ${basedir}/fasta/${name}_macspeaks.dfasta
                	#/bin/mv ${name}_macs_homermotifs_${depth}/background.fa ${basedir}/fasta/${name}_macsbackground.dfasta
                	#gzip ${basedir}/fasta/${name}_macsbackground.dfasta
			#/bin/cp *.bed ${basedir}/bed/
		fi
		#/bin/mv ${name}/*.bedGraph.gz ${basedir}/bedgraph/
		# for histone samples adding option "-size given" (use the whole peak sequence if my guessing is correct) and "-S 15" (for making the search faster) might be useful
		#TODO for naming of the directory instead of using $depth $motiflength would be better, but what if it is like 8,10,13 for example
		# -depth is thrown away in homer2, I am using motiflength(s) for directory naming
	fi
	echo "name=$name iscontrol=$iscontrol"
	if [ $iscontrol ] ; then
		macscontrol="-c ${basedir}/${name}/bam/${name}.bam"
		control="-i ${basedir}/${name}/homer/$name"
		iscontrol=""
	fi
	# --- obsolete, it was for homer1 ---# nice -n 20 findMotifsGenome.pl ${name}_homerpeaks.txt ${genome}r ${name}_homermotifs_$depth -depth allnight -len $motiflength -dumpFasta > ${logdir}/${name}-findMotifsGenome_${depth}.log 2> ${logdir}/${name}-findMotifsGenome_${depth}.err & 
	# allows to run only motif finding after the whole analysis finished
	if [ -f /tmp/${name}bedgraph ] ; then
		/bin/rm /tmp/${name}bedgraph #cleanup, do not leave unused files there
	fi
	echo "finished Homer processing for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
	echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
	# finished HOMER processing, starting ChIPseeqer
	#echo "Starting ChIPseeqer processing" >> ${logdir}/${name}_ChIP-seq_analyze.log
	#cd ${basedir}/${name}/chipseeqer
	#if [ ! -d ${basedir}/${name}/chipseeqer/${name} ] ; then
	#	mkdir ${basedir}/${name}/chipseeqer/${name}
	#	split_bamfile ${basedir}/${name}/bam/${name}.bam -outdir $name > ${logdir}/${name}-split_bamfile.log 2> ${logdir}/${name}-split_bamfile.err 
	#	cd $name
	#	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M ; do 
	#		if [ -f read.chr${i}.bam ] ; then
	#			ln -s read.chr${i}.bam reads.chr${i}
	#		fi
	#	done 
	#	cd ..
	#	ChIPseeqer.bin -chipdir $name -format bam -outfile ${name}_ChIPseeqerpeaks.txt -chrdata $CHIPSEEQERDIR/DATA/${genome}.chrdata > ${logdir}/${name}-ChIPseeqer.log 2> ${logdir}/${name}-ChIPseeqer.err
	#	ChIPseeqerAnnotate --peakfile=${name}_ChIPseeqerpeaks.txt --prefix=${name} --genome=$genome --db=Ensembl > ${logdir}/${name}-ChIPseeqerAnnotate_ensembl.log 2> ${logdir}/${name}-ChIPseeqerAnnotate_ensembl.err
	#	ChIPseeqerAnnotate --peakfile=${name}_ChIPseeqerpeaks.txt --prefix=${name} --genome=$genome --db=RefGene > ${logdir}/${name}-ChIPseeqerAnnotate_refseq.log 2> ${logdir}/${name}-ChIPseeqerAnnotate_refseq.err
	#	ChIPseeqerSummaryPromoters --targets=${name}_ChIPseeqerpeaks.txt --prefix=$name --db=RefGene --genome=$genome --lenu=5000 --lend=1000 >  ${logdir}/${name}-ChIPseeqerSummaryPromoters_refseq.log 2> ${logdir}/${name}-ChIPseeqerSummaryPromoters_refseq.log
	#	ChIPseeqerPeaksTrack --targets=${name}_ChIPseeqerpeaks.txt --trackname="${name}" > ${logdir}/${name}-ChIPseeqerPeaksTrack.log 2> ${logdir}/${name}-ChIPseeqerPeaksTrack.err
	#	ChIPseeqerNongenicAnnotate --peakfile=${name}_ChIPseeqerpeaks.txt --type=RepMasker > ${logdir}/${name}-ChIPseeqerNongenicAnnotate_repeat.log 2> ${logdir}/${name}-ChIPseeqerNongenicAnnotate_repeat.log
	#	ChIPseeqerNongenicAnnotate --peakfile=${name}_ChIPseeqerpeaks.txt --type=CpGislands > ${logdir}/${name}-ChIPseeqerNongenicAnnotate_CpG.log 2> ${logdir}/${name}-ChIPseeqerNongenicAnnotate_CpG.log
	#	ChIPseeqerNongenicAnnotate --peakfile=${name}_ChIPseeqerpeaks.txt --type=SegmentalDups > ${logdir}/${name}-ChIPseeqerNongenicAnnotate_SegmentalDups.log 2> ${logdir}/${name}-ChIPseeqerNongenicAnnotate_SegmentalDups.log
	#	ChIPseeqerFIRE --targets=${name}_ChIPseeqerpeaks.txt --fastafile=${indexdir}/$genome -species=$spec > ${logdir}/${name}-FIRE.log 2> ${logdir}/${name}-FIRE.err
	#	ChIPseeqerMotifMatch --targets=${name}_ChIPseeqerpeaks.txt --suffix=${name}-peaks_vs_$cmotif1 --motifname=${cmotif1}.jaspar -c 2.0 > ${logdir}/${name}-ChIPseeqerMotifMatch_${cmotif1}.log 2>${logdir}/${name}-ChIPseeqerMotifMatch_${cmotif1}.err
	#	ChIPseeqerMotifMatch --targets=${name}_ChIPseeqerpeaks.txt --suffix=${name}-peaks_vs_$cmotif2 --motifname=${cmotif2}.jaspar -c 2.0 > ${logdir}/${name}-ChIPseeqerMotifMatch_${cmotif2}.log 2>${logdir}/${name}-ChIPseeqerMotifMatch_${cmotif2}.err
	#	ChIPseeqeriPAGE --targets=${name}_ChIPseeqerpeaks.txt --lenuP=5000 --lendP=1000 --suffix=${name}_PAGE --genome=$genome --db=RefGene > ${logdir}/${name}-PAGE.log 2> ${logdir}/${name}-PAGE.err
	#	if [ -d  $phastdir -a "$species" = "hs" ] ; then
	#		ChIPseeqerCons --targets=${name}_ChIPseeqerpeaks.txt --consdir $phastdir --outfile=${name}_TF_targets_conservation.txt --outepsmap=${name}_TF_targets_conservation.eps --format=gzscores --category=placental --make_rand=1 --show_profiles=1 --around_summit=1 --showalldata=1
	#		ChIPseeqerFindDistalPeaks --targets=${name}_ChIPseeqerpeaks.txt --db=RefGene --mindistaway=10000 --prefix=${name}_distpeaks --genome=$genome >${logdir}/${name}-ChIPseeqerFindDistalPeaks.log 2> ${logdir}/${name}-ChIPseeqerFindDistalPeaks.err > ${logdir}/${name}-ChIPseeqerCons.log 2> ${logdir}/${name}-ChIPseeqerCons.err
	#		$CHIPSEEQERDIR/SCRIPTS/extract_regions_around_peak_summits.pl --peakfile=${name}_distpeaks.RefGene.DISTPEAKS --w=2000 > ${logdir}/${name}-extract_regions_around_peak_summits.log 2> ${logdir}/${name}-extract_regions_around_peak_summits.err
	#		ChIPseeqer2Cons.bin -regions ${name}_distpeaks.RefGene.DISTPEAKS.2kbaround --consdir $phastdir  -format gzscores -category placental \
	#		-method phastCons -make_rand 1 -randist 50000 -outrandom ${name}_distpeaks.RefGene.DISTPEAKS.2kbaround.randcons -show_profiles 1 \
	#		-outfile ${name}_distpeaks_out.cons -outepsmap ${name}_distpeaks_out.eps > ${logdir}/${name}-ChIPseeqer2Cons.log 2> ${logdir}/${name}-ChIPseeqer2Cons.err
	#	fi
	#	/bin/cp *.wgl.gz ${basedir}/bed/
	#	echo "finished ChIPseeqer processing for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
	#else
	#	echo "ChIPseeqer processing for $name already done" >> ${logdir}/${name}_ChIP-seq_analyze.log
	#	# Put here other ChIPseeqer programs to run or re-run after the basic analysis has been finished
	#	ChIPseeqeriPAGE --targets=${name}_ChIPseeqerpeaks.txt --lenuP=5000 --lendP=1000 --suffix=${name}_PAGE --genome=$genome --db=RefGene > ${logdir}/${name}-PAGE.log 2> ${logdir}/${name}-PAGE.err
	#fi
	echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
	# finished ChIPseeqer processing
	#echo "Starting MEME de novo motif finding" >> ${logdir}/${name}_ChIP-seq_analyze.log
	#cd ${basedir}/${name}/meme
	# meme paralell or meme for one processor. uncomment the appropiate one to run (it can be very slow)
	#mpirun -np $procno meme_p ${basedir}/${name}/homer/${name}_homermotifs_$depth/target.fa -dna -mod zoops -revcomp -oc ${minmotifl}_${maxmotifl}_p -nmotifs $nmotif -minw $minmotifl -maxw $maxmotifl < /dev/null > ${logdir}/${name}-meme_${minmotifl}_${maxmotifl}_p.log 2> ${logdir}/${name}-meme_${minmotifl}_${maxmotifl}_p.err
	#meme_p ${basedir}/${name}/homer/${name}_homermotifs_$depth/target.fa -dna -mod zoops -revcomp -oc ${minmotifl}_$maxmotifl -nmotifs $nmotif -minw $minmotifl -maxw $maxmotifl < /dev/null > ${logdir}/${name}-meme_${minmotifl}_${maxmotifl}.log 2> ${logdir}/${name}-meme_${minmotifl}_${maxmotifl}.err
	#echo "finished MEME processing for $name at `date`" >> ${logdir}/${name}_ChIP-seq_analyze.log
	#echo "------------------------------------------------------------" >> ${logdir}/${name}_ChIP-seq_analyze.log
done
