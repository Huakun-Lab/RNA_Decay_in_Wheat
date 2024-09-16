###########################################
#                                         #
#Pearson correlation coefficient          #
#Bivariate Correlations analysis          #
###########################################

library(ggstatsplot)
options(stringsAsFactors=F)


DR1 <-read.table("input.txt",header = T,sep="\t")
cor.test(DR1$variables1,DR1$variables2, method = "pearson")


pdf("output2.pdf", width = 5.0, height = 3.2)
ggscatterstats(DR1, x = "variables1",
               y = "variables2", 
               title = "Title")
dev.off()
