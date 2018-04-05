require(grid)
require(ggplot2)
setwd("/data6/working_groups/nlab_group/IL4_STAT6/ATAC-Seq/closed_and_open/closed_PU1_project/50M_PU1_v9/TF_map")

tt <- read.table("all_vs_alone.tsv"); tt <- cbind(Name = rownames(tt), tt); tt <- cbind(tt, Percent = tt$Alone/tt$All)
#arr <- arrow(ends = ifelse(factor2.y<=factor1.y, "first", "last"), length = unit(0.05, "inches")); null.arr <- NULL
#curvs <- matrix(c(NA, NA, NA, NA, NA,
#		  0.24, NA,  NA, NA,  NA,
#		  0.725, -0.23,   NA,  NA,  NA,
#		  0.63, 1.35,   -0.25,  NA,   NA,
#		  0.73, 1.25,   0.35,  0.65,   NA
#), nrow = 5, byrow = T)
curvs <- matrix(c(NA, 0.2, 0.7, 0.6, 0.5,
                  NA, NA,   -0.35, -0.5,  1.15,
                  NA, NA,   NA,  -0.28,  0.43,
                  NA, NA,   NA,  NA,   0.5,
                  NA, NA,   NA,  NA,   NA
), nrow = 5, byrow = T)

overlaps <- as.matrix(read.table("percents.tsv"))
t0 <- as.data.frame(matrix(c(0,0,0,65000,65000,1), byrow = T, nrow = 2, dimnames = list(c("x", "y"), c("All", "Alone", "Percent"))))
#for (i in 1:5) 
#	percents[i,] <- overlaps[i,]/diag(overlaps)[i]
colnames(curvs) <- rownames(curvs) <- c("PU1", "IRF8", "JUNB", "RUNX1", "CEBPa")
gg <- ggplot(t0, aes(Alone, All)) + geom_line(data = t0, linetype = 2, alpha = 0)
for (i in 1:ncol(tt)) {
#for (i in 1) {
	for (j in 1:nrow(tt)) {
		factor1 <- tt$Name[i]; factor2 <- tt$Name[j]
		if (factor1 == factor2 || is.na(curvs[i,j])) next;
		factor1.x <- tt[factor1==tt[,1],"Alone"]; factor1.y <- tt[factor1==tt[,1],"All"]
		factor2.x <- tt[factor2==tt[,1],"Alone"]; factor2.y <- tt[factor2==tt[,1],"All"]
		pp <- overlaps[i,j]/35000
		lwd <- pp*20
		col <- ifelse(factor2.y<=factor1.y, "blue", "blue")
		cat("factor1=", tt$Name[i], factor1.x, factor1.y, "factor2=", tt$Name[j], factor2.x, factor2.y, "curv=", curvs[i,j],"pp=",pp, "i=", i, "j=",j,"\n") 
		curve1 <- curveGrob(0, 0, 1, 1, default.units = "npc", curvature = curvs[i,j], angle = 90, ncp = 20, shape = 1, square = F, squareShape = 1,inflect = F, arrow = null.arr, open = T, debug = F, name = NULL, gp = gpar(col = col, lwd = lwd), vp = NULL)
		print(lwd)	
		gg <- gg + annotation_custom(grob = curve1, factor1.x, factor2.x, ifactor1.y, factor2.y) 
	}
}
gg.Percent <- gg + geom_point(data = tt, aes(Alone, All), size = tt$All/1500, alpha = 1, color = c("#FF0000", rep("#006666", 4)), shape = 20)+ xlim(0,45000) + ylim(0,70000) + theme_bw()
ggsave(gg.Percent, filename = "bluess_abs_final.pdf", width = 8, height = 8)
