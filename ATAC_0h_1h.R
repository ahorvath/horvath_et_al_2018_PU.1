library(DiffBind)
options(scipen = 999)
# create result dir (if it is necessary)
result.dir <- "results_FC1_pv0.05"

# set work dir
if (!file.exists(result.dir)) {
	dir.create(result.dir)
}

setwd(result.dir)

# read raw data
ATAC_0h_1h.dba = dba(sampleSheet = "csv/ATAC_0h_1h.csv", peakCaller = "bed")

# create consensus peak set
ATAC_0h_1h.dba.count = dba.count(ATAC_0h_1h.dba)

print(ATAC_0h_1h.dba)
print(ATAC_0h_1h.dba.count)

# draw heatmap plot based on full raw data
#pdf("ATAC_0h_1h_full_heatmap.pdf")
#dba.plotHeatmap(ATAC_0h_1h.dba, correlations = F, margin = 15)
#dev.off()

#pdf("ATAC_0h_1h_full_heatmap_consensus_peak.pdf")
#dba.plotHeatmap(ATAC_0h_1h.dba.count, correlations = F, margin = 15, maxval = 100)
#dev.off()

# draw and save correlation plot based on raw data
#pdf("ATAC_0h_1h_dba_plot_before.pdf")
#plot(ATAC_0h_1h.dba, margin = 15)
#dev.off()

# draw and save correlation plot based on consensus peak set
#pdf("ATAC_0h_1h_dba_plot_after.pdf")
#plot(ATAC_0h_1h.dba.count, margin = 15, colScheme="Blues")
#dev.off()

# create contrast(s)
ATAC_0h_1h.contrast = dba.contrast(ATAC_0h_1h.dba.count, categories = DBA_CONDITION, minMembers = 2)
ATAC_0h_1h.analyze = dba.analyze(ATAC_0h_1h.contrast, bFullLibrarySize = T) #, DBA_DESEQ,DBA_EDGER_CLASSIC,DBA_DESEQ_CLASSIC,DBA_EDGER_GLM,DBA_DESEQ_GLM))
th <- 0.05
# create correlation heatmap based on changing peak (p-value)
#pdf("ATAC_0h_1h_corr_heatmap.pdf")
#dba.plotHeatmap(ATAC_0h_1h.analyze, correlation = T, th = th, bUsePval = T, contrast = 1, margin = 15,  colScheme="Blues")
#dev.off()

# create heatmap based on changing peak (simple p-value)
#pdf("ATAC_0h_1h_heatmap.pdf")
#dba.plotHeatmap(ATAC_0h_1h.analyze, correlation = F, th = th, bUsePval = T, contrast = 1, margin = 15,  colScheme="Blues")
#dev.off()
FC <- 1
# create expression table
ATAC_0h_1h.DB = dba.report(ATAC_0h_1h.analyze, th = th, contrast = 1, bUsePval = T, fold = log2(FC))
ATAC_0h_1h.table <- as.data.frame(ATAC_0h_1h.DB)
ATAC_0h_1h.table <- ATAC_0h_1h.table[order(ATAC_0h_1h.table[,1]),]

ATAC_0h_1h.table.up <- ATAC_0h_1h.table[ATAC_0h_1h.table$Fold < 0,]
ATAC_0h_1h.table.up.poised <- subset(ATAC_0h_1h.table.up, Conc_0h <= 1)
ATAC_0h_1h.table.up.const <- subset(ATAC_0h_1h.table.up, Conc_0h > 1)
ATAC_0h_1h.table.down <- ATAC_0h_1h.table[ATAC_0h_1h.table$Fold > 0,]

# create binding affinity boxplot based on changing peak (simple p-value)

#dev.off()

# write expression tables
write.table(cbind(ATAC_0h_1h.table.down[,1:3], paste(paste(ATAC_0h_1h.table.down[,1], ATAC_0h_1h.table.down[,2], sep = ":"), ATAC_0h_1h.table.down[,3], sep = "-"), ATAC_0h_1h.table.down[,9]),"ATAC_0h_1h_pv_0.05_FC1_down_regulated.bed", row.names = F, col.names = F, quote = F, sep = "\t")

head(ATAC_0h_1h.table.up)
write.table(cbind(ATAC_0h_1h.table.up[,1:3], paste(paste(ATAC_0h_1h.table.up[,1], ATAC_0h_1h.table.up[,2], sep = ":"), ATAC_0h_1h.table.up[,3], sep = "-"), ATAC_0h_1h.table.up[,9]),"ATAC_0h_1h_pv_0.05_FC1_up_regulated.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(cbind(ATAC_0h_1h.table.up.poised, paste(paste(ATAC_0h_1h.table.up.poised[,1], ATAC_0h_1h.table.up.poised[,2], sep = ":"), ATAC_0h_1h.table.up.poised[,3], sep = "-"), ATAC_0h_1h.table.up.poised[,9]),"ATAC_0h_1h_up_regulated_poised.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(cbind(ATAC_0h_1h.table.up.const, paste(paste(ATAC_0h_1h.table.up.const[,1], ATAC_0h_1h.table.up.const[,2], sep = ":"), ATAC_0h_1h.table.up.const[,3], sep = "-"), ATAC_0h_1h.table.up.const[,9]),"ATAC_0h_1h_up_regulated_const.bed", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(cbind(ATAC_0h_1h.table, paste(paste(ATAC_0h_1h.table[,1], ATAC_0h_1h.table[,2], sep = ":"), ATAC_0h_1h.table[,3], sep = "-"), ATAC_0h_1h.table[,9]),"diff_bind_ATAC_0h_1h_IL4_1h_sig.bed", row.names = F, col.names = F, quote = F, sep = "\t")

ATAC_0h_1h.DB = dba.report(ATAC_0h_1h.analyze, th = 1, contrast = 1, bUsePval = T, fold = log2(1))
ATAC_0h_1h.table <- as.data.frame(ATAC_0h_1h.DB)
ATAC_0h_1h.full.table <- ATAC_0h_1h.table[order(ATAC_0h_1h.table[,1]),]

write.table(ATAC_0h_1h.full.table, "diff_bind_ATAC_0h_1h_full_table.bed", row.names = F, col.names = F, quote = F, sep = "\t")
save.image()
