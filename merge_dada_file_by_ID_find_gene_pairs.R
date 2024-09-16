###############################################
#                                             #
#The task is to find the rows with the        #
#same gene IDs in two input files, merge      #
#them, and output the merged result. Then,    #
#based on the homoeologous gene set of the    #
#AB part, perform homoeologous gene pairing.  #
###############################################

library(dplyr)
options(stringsAsFactors=F)

df1 <- read.delim("input1.txt",header = T,sep="\t")

df2 <- read.delim("input2.txt",header = T,sep="\t")

df1 %>% dplyr::mutate(rate_log=(log10(alpha)))-> df1

df3 <- unique(merge(df1,df2,by="gene_id"))
df3 <- filter(df3,df3["group"]=="all")
colnames(df3)[1]<- "gene_id_A"

df4 <- read.table("AvsB_gene_pairs.txt",header = T,sep="\t")
df5 <- unique(merge(df3,df4,by="gene_id_A"))
colnames(df3)[1]<- "gene_id_B"
df6 <- unique(merge(df5,df3,by="gene_id_B"))

write.table(df6,file="gene_pairs_containg_rate_rpm_info.txt", quote = FALSE,row.names = T, sep = "\t")


