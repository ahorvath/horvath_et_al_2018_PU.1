library(ggplot2)
 
# Create test data.
dat = data.frame(count=c(29539, 3566, 1879), category=c("NC", "Down", "Up"))
 
# Add addition columns, needed for drawing with geom_rect.
dat$fraction = dat$count / sum(dat$count)
dat = dat[order(dat$fraction), ]
dat$ymax = cumsum(dat$fraction)
dat$ymin = c(0, head(dat$ymax, n=-1))
dat$category <- factor(dat$category, levels = c("NC", "Down", "Up"))
# Make the plot
p1 = ggplot(dat, aes(fill=category, ymax=ymax, ymin=ymin, xmax=3, xmin=2)) +
     geom_rect() +
     coord_polar(theta="y") +
     xlim(c(0, 3)) +
     theme(panel.grid=element_blank()) +
     theme(axis.text=element_blank()) +
     theme(axis.ticks=element_blank()) +
     theme_bw() +
     annotate("text", x = 0, y = 0, label = "My Ring plot !") +
     labs(title="") + 
     scale_fill_manual(values = c("grey", "blue", "red"))
p1

ggsave(p1, filename = "donuts.jpg")
