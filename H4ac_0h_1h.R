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
H4ac_0h_1h.dba = dba(sampleSheet = "/data6/working_groups/nlab_group/IL4_STAT6/ATAC-Seq/closed_and_open/closed_PU1_project/50M_PU1_v7/DiffBind/csv/H4ac_0h_1h.csv", caller = "bed")

# create consensus peak set
H4ac_0h_1h.dba.count = dba.count(H4ac_0h_1h.dba)
print(H4ac_0h_1h.dba)
print(H4ac_0h_1h.dba.count)

# draw heatmap plot based on full raw data
pdf("H4ac_0h_1h_full_heatmap.pdf")
dba.plotHeatmap(H4ac_0h_1h.dba, correlations = F, margin = 15)
dev.off()

pdf("H4ac_0h_1h_full_heatmap_consensus_peak.pdf")
dba.plotHeatmap(H4ac_0h_1h.dba.count, correlations = F, margin = 15, maxval = 100)
dev.off()

# draw and save correlation plot based on raw data
pdf("H4ac_0h_1h_dba_plot_before.pdf")
plot(H4ac_0h_1h.dba, margin = 15,  colScheme="Purples")
dev.off()

# draw and save correlation plot based on consensus peak set
pdf("H4ac_0h_1h_dba_plot_after.pdf")
plot(H4ac_0h_1h.dba.count, margin = 15,  colScheme="Purples")
dev.off()

# create contrast(s)
H4ac_0h_1h.contrast = dba.contrast(H4ac_0h_1h.dba.count, categories = DBA_CONDITION, minMembers = 2)
H4ac_0h_1h.analyze = dba.analyze(H4ac_0h_1h.contrast, bFullLibrarySize = T) #, DBA_DESEQ,DBA_EDGER_CLASSIC,DBA_DESEQ_CLASSIC,DBA_EDGER_GLM,DBA_DESEQ_GLM))
th <- 0.05
# create correlation heatmap based on changing peak (p-value)
pdf("H4ac_0h_1h_corr_heatmap.pdf")
dba.plotHeatmap(H4ac_0h_1h.analyze, correlation = T, th = th, bUsePval = T, contrast = 1, margin = 15,  colScheme="Purples")
dev.off()

# create heatmap based on changing peak (simple p-value)
pdf("H4ac_0h_1h_heatmap.pdf")
dba.plotHeatmap(H4ac_0h_1h.analyze, correlation = F, th = th, bUsePval = T, contrast = 1, margin = 15,  colScheme="Purples")
dev.off()

# create expression table
H4ac_0h_1h.DB = dba.report(H4ac_0h_1h.analyze, th = th, contrast = 1, bUsePval = T, fold = log2(1))
H4ac_0h_1h.table <- as.data.frame(H4ac_0h_1h.DB)
H4ac_0h_1h.table <- H4ac_0h_1h.table[order(H4ac_0h_1h.table[,1]),]

H4ac_0h_1h.table.up <- H4ac_0h_1h.table[H4ac_0h_1h.table$Fold < 0,]
H4ac_0h_1h.table.up.poised <- subset(H4ac_0h_1h.table.up, Conc_0h <= 7)
H4ac_0h_1h.table.up.const <- subset(H4ac_0h_1h.table.up, Conc_0h > 7)
H4ac_0h_1h.table.down <- H4ac_0h_1h.table[H4ac_0h_1h.table$Fold > 0,]

# create binding affinity boxplot based on changing peak (simple p-value)
pdf("H4ac_0h_1h_dba_plot_box.pdf")
dba.plotBox(H4ac_0h_1h.analyze, th = th, bUsePval = T, contrast = 1, margin = 15, fold = log2(1))
dev.off()

# write expression tables
write.table(cbind(H4ac_0h_1h.table.down[,c(1,2,3)], paste(paste(H4ac_0h_1h.table.down[,1], H4ac_0h_1h.table.down[,2], sep = ":"), H4ac_0h_1h.table.down[,3], sep = "-"), H4ac_0h_1h.table.down[,9]),"H4ac_0h_1h_down_regulated.bed", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(cbind(H4ac_0h_1h.table.up[,c(1,2,3)], paste(paste(H4ac_0h_1h.table.up[,1], H4ac_0h_1h.table.up[,2], sep = ":"), H4ac_0h_1h.table.up[,3], sep = "-"), H4ac_0h_1h.table.up[,9]),"H4ac_0h_1h_up_regulated.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(cbind(H4ac_0h_1h.table.up.poised[,c(1,2,3)], paste(paste(H4ac_0h_1h.table.up.poised[,1], H4ac_0h_1h.table.up.poised[,2], sep = ":"), H4ac_0h_1h.table.up.poised[,3], sep = "-"), H4ac_0h_1h.table.up.poised[,9]),"H4ac_0h_1h_up_regulated_poised.bed", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(cbind(H4ac_0h_1h.table.up.const[,c(1,2,3)], paste(paste(H4ac_0h_1h.table.up.const[,1], H4ac_0h_1h.table.up.const[,2], sep = ":"), H4ac_0h_1h.table.up.const[,3], sep = "-"), H4ac_0h_1h.table.up.const[,9]),"H4ac_0h_1h_up_regulated_const.bed", row.names = F, col.names = F, quote = F, sep = "\t")

write.table(cbind(H4ac_0h_1h.table[,c(1,2,3)], paste(paste(H4ac_0h_1h.table[,1], H4ac_0h_1h.table[,2], sep = ":"), H4ac_0h_1h.table[,3], sep = "-"), H4ac_0h_1h.table[,9]),"diff_bind_H4ac_0h_1h_IL4_1h_sig.bed", row.names = F, col.names = F, quote = F, sep = "\t")

H4ac_0h_1h.full.table <- as.data.frame(H4ac_0h_1h.dba.count$peaks[1])[,c(1:3)]
write.table(cbind(H4ac_0h_1h.full.table[, c(1:3)], paste(paste(H4ac_0h_1h.full.table[, 1],H4ac_0h_1h.full.table[, 2], sep = ":"), H4ac_0h_1h.full.table[, 3], sep = "-")), "diff_bind_H4ac_0h_1h_full_table.bed", row.names = F, col.names = F, quote = F, sep = "\t")
save.image()
