#/usr/local/src/R-2.14.1/bin/R
#library(DiffBind,lib.loc="/molbio/flatfiles/R-alternative_library/DiffBind-1.2.4/")
library(DiffBind)
setwd("/data6/working_groups/nlab_group/IL4_STAT6/ATAC-Seq/closed_and_open/closed_PU1_project/50M_PU1_v9/NFR_vs_Factor_count/DiffBind")
min.members <- 1:5
for (min.member in min.members) {
	print(min.member)
#	PU1.dba = dba(sampleSheet = "csv/diffbind_all.csv", caller = "bed", minOverlap = min.member)
	PU1.dba = dba(sampleSheet = "csv/diffbind_all.csv", peakFormat = "bed", minOverlap = min.member)
#	PU1.dba = dba(sampleSheet = "csv/diffbind_all.csv", minOverlap = min.member, peakFormat = "bed", peakCaller = "macs", scoreCol = 5)
	PU1.dba.count = dba.count(PU1.dba, bLog = T, minOverlap = min.member)
	print(PU1.dba.count)
	full.PU1.table <- as.data.frame(PU1.dba.count$peaks)[,1:3]

	write.table(full.PU1.table,paste("consensus_min", min.member, ".bed", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
}
