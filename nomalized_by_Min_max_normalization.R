###########################################
#                                         #
#Min-max normalization                    #
#                                         #
###########################################


library(dplyr)

RPMs_mean <- read.delim("RPMs_mean_2_replicates_relative_T0.txt",header = T,sep="\t") 
RPMs_mean2 <-  as.matrix(RPMs_mean[,2:9])


min.max.norm <- function(x){
  ((x-min(x))/(max(x)-min(x)))
}

RPMs_nor <- apply(RPMs_mean2[1:43670,],1,min.max.norm)
RPMs_nor2 <- round(RPMs_nor,4)

head(RPMs_nor2)
class(RPMs_nor2)

RPMs_nor3 <- t(RPMs_nor2)
write.table(RPMs_nor3,file="RPMs_relative_nomalized_by_0_1_output.txt", quote = FALSE,row.names = F, sep = "\t")














