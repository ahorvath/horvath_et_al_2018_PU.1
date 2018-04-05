library(reshape2) 
library(ggplot2)
library(RColorBrewer)
library(limma)

#margin <- "_distal_summits_pm100"
margin <- "_distal"
	tt <- NULL
	ttt <- NULL
	for (file in list.files(full.names = F,  pattern = paste0(".*_on_w_old_up.*_distal.tsv"))) {
		print(margin)
		current.tt <- read.table(file, header = F, sep = "\t")
		tt <- rbind(tt, cbind(current.tt[,4, drop = F], Peakset = gsub(".tsv", "", file)))
		ttt <- cbind(ttt, current.tt[,4])
	}
	colnames(ttt) <- gsub("mm_BMDM_|_distal.tsv", "", list.files(full.names = F,  pattern = paste0(".*_on_w_old_down.*_distal.tsv")))

	melted.tt <- data.frame(matrix(unlist(strsplit(x = as.character(tt$Peakset), split = "_on_")), ncol = 2, byrow=T))

	melted.tt <- cbind(melted.tt, tt$V4)
	colnames(melted.tt) <- c("Sample", "Peakset", "Occupancy")
	melted.tt$Peakset <- as.factor(gsub("_distal", "", melted.tt$Peakset))
#	melted.tt$Peakset <- gsub(paste0(margin, "|", name, "_on_"), "", melted.tt$Peakset)
#	melted.tt$Peakset <-  factor(melted.tt$Peakset, levels = c("PU1_0h_not_ATAC_0h", "PU1_0h_and_ATAC_0h", "ATAC_0h_not_PU1_0h"))

	gg <- ggplot(melted.tt, aes(Peakset, log10(Occupancy), fill = Sample)) + geom_boxplot(alpha = 1, outlier.shape = 1) + scale_fill_manual(values = c("#FF0000", "#0000FF")) + theme(legend.title = element_blank(), axis.title = element_text(size = 20), legend.text = element_text(size = 20), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() + ylim(-1.0, 2.5)
	name <- "w_old_up_boxplots"
	ggsave(gg, filename = paste0(name,"_Occupancies.pdf"), width = 5, height = 6)

	gg <- ggplot(melted.tt, aes(Peakset, log10(Occupancy), fill = Peakset)) + geom_boxplot(alpha = 1, outlier.shape = 1) + scale_fill_manual(values = rep("#0000FF",3)) + theme(legend.title = element_blank(), axis.title = element_text(size = 20), legend.text = element_text(size = 20), strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), axis.text.x = element_text(angle = 90, hjust = 1)) + theme_bw() + ylim(-1.0, 2.5) + theme(axis.title=element_blank(), axis.text=element_blank(), legend.position="none")
	ggsave(gg, filename = paste0(name, margin, "_Occupancies_no_labels.pdf"), width = 2.5, height = 6)
