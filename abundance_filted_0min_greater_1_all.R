###########################################
#                                         #
#Low-abundance transcripts were filted    # 
#based on three biological replicates     #
#of 0min data.                            #
###########################################

library(dplyr)
options(stringsAsFactors=F)


df1 <- read.table("inputfile.txt",header = T,sep="\t") 

df2 <- unique(filter(df1,df1$kro_00_r1> 1 & df1$kro_00_r2 > 1 & df1$kro_00_r3 > 1))

write.table(df2, file="output_filted_file.txt", quote = FALSE,row.names = T,sep = "\t",fileEncoding="UTF-8")


