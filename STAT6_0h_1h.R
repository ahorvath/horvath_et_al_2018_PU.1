#/usr/local/src/R-2.14.1/bin/R
library(DiffBind,lib.loc="/molbio/flatfiles/R-alternative_library/DiffBind-1.2.4/")

# create result dir (if it is necessary)
result.dir <- "/data6/working_groups/nlab_group/IL4_STAT6/ATAC-Seq/closed_and_open/closed_PU1_project/50M_PU1_v7/DiffBind/results"

# set work dir
if (!file.exists(result.dir)) {
	dir.create(result.dir)
}

setwd(result.dir)

# read raw data
STAT6_0h_1h.dba = dba(sampleSheet = "/data6/working_groups/nlab_group/IL4_STAT6/ATAC-Seq/closed_and_open/closed_PU1_project/50M_PU1_v7/DiffBind/csv/STAT6_0h_1h.csv", caller = "bed")

# create consensus peak set
STAT6_0h_1h.dba.count = dba.count(STAT6_0h_1h.dba)

print(STAT6_0h_1h.dba)
print(STAT6_0h_1h.dba.count)

# draw heatmap plot based on full raw data
pdf("STAT6_0h_1h_full_heatmap.pdf")
dba.plotHeatmap(STAT6_0h_1h.dba, correlations = F, margin = 15)
dev.off()

pdf("STAT6_0h_1h_full_heatmap_consensus_peak.pdf")
dba.plotHeatmap(STAT6_0h_1h.dba.count, correlations = F, margin = 15, maxval = 100)
dev.off()

# draw and save correlation plot based on raw data
pdf("STAT6_0h_1h_dba_plot_before.pdf")
plot(STAT6_0h_1h.dba, margin = 15)
dev.off()

# draw and save correlation plot based on consensus peak set
pdf("STAT6_0h_1h_dba_plot_after.pdf")
plot(STAT6_0h_1h.dba.count, margin = 15)
dev.off()

# create contrast(s)
STAT6_0h_1h.contrast = dba.contrast(STAT6_0h_1h.dba.count, categories = DBA_CONDITION, minMembers = 2)
STAT6_0h_1h.analyze = dba.analyze(STAT6_0h_1h.contrast, bFullLibrarySize = T) #, DBA_DESEQ,DBA_EDGER_CLASSIC,DBA_DESEQ_CLASSIC,DBA_EDGER_GLM,DBA_DESEQ_GLM))
th <- 0.05

# create correlation heatmap based on changing peak (p-value)
pdf("STAT6_0h_1h_corr_heatmap.pdf")
dba.plotHeatmap(STAT6_0h_1h.analyze, correlation = T, th = th, bUsePval = T, contrast = 1, margin = 15,  colScheme="Greens")
dev.off()

# create heatmap based on changing peak (simple p-value)
pdf("STAT6_0h_1h_heatmap.pdf")
dba.plotHeatmap(STAT6_0h_1h.analyze, correlation = F, th = th, bUsePval = T, contrast = 1, margin = 15,  colScheme="Greens")
dev.off()
FC <- 1.5
# create expression table
STAT6_0h_1h.DB = dba.report(STAT6_0h_1h.analyze, th = th, contrast = 1, bUsePval = T, fold = log2(FC))
STAT6_0h_1h.table <- as.data.frame(STAT6_0h_1h.DB)
STAT6_0h_1h.table <- STAT6_0h_1h.table[order(STAT6_0h_1h.table[,1]),]

STAT6_0h_1h.table.up <- STAT6_0h_1h.table[STAT6_0h_1h.table$Fold < 0,]
STAT6_0h_1h.table.down <- STAT6_0h_1h.table[STAT6_0h_1h.table$Fold > 0,]
# create binding affinity boxplot based on changing peak (simple p-value)
pdf("STAT6_0h_1h_dba_plot_box.pdf")
dba.plotBox(STAT6_0h_1h.analyze, th = th, bUsePval = T, contrast = 1, margin = 15, fold = log2(FC))
dev.off()


# write expression tables
write.table(cbind(STAT6_0h_1h.table.down[,c(1,2,3)], paste(paste(STAT6_0h_1h.table.down[,1], STAT6_0h_1h.table.down[,2], sep = ":"), STAT6_0h_1h.table.down[,3], sep = "-"), STAT6_0h_1h.table.down[,9]),"STAT6_0h_1h_down_regulated.bed", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(cbind(STAT6_0h_1h.table.up[,c(1,2,3)], paste(paste(STAT6_0h_1h.table.up[,1], STAT6_0h_1h.table.up[,2], sep = ":"), STAT6_0h_1h.table.up[,3], sep = "-"), STAT6_0h_1h.table.up[,9]),"STAT6_0h_1h_up_regulated.bed", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(cbind(STAT6_0h_1h.table[,c(1,2,3)], paste(paste(STAT6_0h_1h.table[,1], STAT6_0h_1h.table[,2], sep = ":"), STAT6_0h_1h.table[,3], sep = "-"), STAT6_0h_1h.table[,9]),"diff_bind_STAT6_0h_1h_IL4_1h_sig.bed", row.names = F, col.names = F, quote = F, sep = "\t")

STAT6_0h_1h.full.table <- as.data.frame(STAT6_0h_1h.dba.count$peaks[1])[,c(1:3)]
write.table(cbind(STAT6_0h_1h.full.table[, c(1:3)], paste(paste(STAT6_0h_1h.full.table[, 1],STAT6_0h_1h.full.table[, 2], sep = ":"), STAT6_0h_1h.full.table[, 3], sep = "-")), "diff_bind_STAT6_0h_1h_full_table.bed", row.names = F, col.names = F, quote = F, sep = "\t")
save.image()
