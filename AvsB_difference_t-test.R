###########################################
#                                         #
#Pearson correlation coefficient          #
#Bivariate Correlations analysis          #
###########################################

library("dplyr")

DEN<- read.delim("Rate_pairs_input1.txt",header = T,sep="\t")


DEN%>%dplyr::mutate(p-value=apply(DEN,1,function(x) t.test(as.numeric(x[6:9]),as.numeric(x[10:13]),paired = TRUE)$p.value))->DEN


write.table(DEN,file = "Rate_pairs_t_test_output.txt",row.names = FALSE,sep="\t")
