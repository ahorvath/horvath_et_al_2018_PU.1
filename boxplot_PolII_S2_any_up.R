library(ggplot2)

numbers <- read.table("../numbers.tsv")


denovo.veh <- read.table("boxplots/mm_BMDM_vehicle_PoI_II_S2_CS1361_S1_on_denovo_PolII_any_up_pos.tsv")
labeled.veh <- read.table("boxplots/mm_BMDM_vehicle_PoI_II_S2_CS1361_S1_on_labeled_PolII_any_up_pos.tsv")
poised.veh <- read.table("boxplots/mm_BMDM_vehicle_PoI_II_S2_CS1361_S1_on_poised_PolII_any_up_pos.tsv")
const.veh <- read.table("boxplots/mm_BMDM_vehicle_PoI_II_S2_CS1361_S1_on_const_PolII_any_up_pos.tsv")

denovo.IL4 <- read.table("boxplots/mm_BMDM_IL4_PoI_II_S2_CS1362_S1_on_denovo_PolII_any_up_pos.tsv")
labeled.IL4 <- read.table("boxplots/mm_BMDM_IL4_PoI_II_S2_CS1362_S1_on_labeled_PolII_any_up_pos.tsv")
poised.IL4 <- read.table("boxplots/mm_BMDM_IL4_PoI_II_S2_CS1362_S1_on_poised_PolII_any_up_pos.tsv")
const.IL4 <- read.table("boxplots/mm_BMDM_IL4_PoI_II_S2_CS1362_S1_on_const_PolII_any_up_pos.tsv")

tt <- as.data.frame(
	rbind(cbind(denovo.veh[,4], Type = "denovo", Treat = "veh"),
	     cbind(labeled.veh[,4], Type = "labeled", Treat = "veh"),
	     cbind(poised.veh[,4], Type = "poised", Treat = "veh"),
	     cbind(const.veh[,4], Type = "const", Treat = "veh"),
	     cbind(denovo.IL4[,4], Type = "denovo", Treat = "IL4_1h"),
	     cbind(labeled.IL4[,4], Type = "labeled", Treat = "IL4_1h"),
	     cbind(poised.IL4[,4], Type = "poised", Treat = "IL4_1h"),
	     cbind(const.IL4[,4], Type = "const", Treat = "IL4_1h")
))


denovo.fold <- mean((denovo.IL4[,4]+0.1)/(denovo.veh[,4]+0.1))
labeled.fold <- mean((labeled.IL4[,4]+0.1)/(labeled.veh[,4]+0.1))
poised.fold <- mean((poised.IL4[,4]+0.1)/(poised.veh[,4]+0.1))
const.fold <- mean((const.IL4[,4]+0.1)/(const.veh[,4]+0.1))

colnames(tt)[1] <- "log10_RPKM"
tt[,1] <- log10(as.numeric(as.character(tt[,1])))
tt$Treat <- factor(tt$Treat, levels = c("veh", "IL4_1h"))
tt$Type <- factor(tt$Type, levels = c("denovo", "labeled", "poised", "const"))
gg <- ggplot(tt, aes(Treat, log10_RPKM, fill = Treat)) + geom_boxplot() + ylim(-0.5,1) + facet_wrap( ~ Type, ncol = 4) + theme_bw() +  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + scale_fill_manual(values = c("#F2F0F7", "#756BB1"))
ggsave(gg, filename = "boxplot_PolII_S2_on_any_subclasses.pdf", width = 8)

