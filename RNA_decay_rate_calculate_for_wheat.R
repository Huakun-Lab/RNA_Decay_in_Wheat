##########################################################
#Maximum Likelihood Decay Modeling of RNA Decay Data     #
#                                                        #
#ML model test method                                    #
#BiocManager::install("RNAdecay")               .        #
#                                                        #
#reference:                                              #
#Sorenson R S, Deshotel M J, Johnson K, et al.PNAS,      #
#2018,                                                   #
# https://doi.org/10.1073/pnas.1712312115.               #
# https://doi.org/doi:10.18129/B9.bioc.RNAdecay          #
##########################################################

library(RNAdecay)
library(dplyr)

#K_RPM <- read.delim("kro_abundance_add_WT.txt",header = T,sep = "\t", fileEncoding="UTF-8")
#K_RPM_1<- filter(K_RPM,(K_RPM$kro_00_r1 >1 &  K_RPM$kro_00_r2 >1 & K_RPM$kro_00_r3 >1 ) )
#write.table(K_RPM_1,file="input.txt", quote = FALSE,row.names = T,sep = "\t",fileEncoding="UTF-8")

RPMs <- read.table("input.txt",header = T) 

cols(patterns = c("WT", "00"), df = RPMs) 

colnames(RPMs)[cols(patterns = c("WT", "00"), RPMs, "WT", "00")]

wrdir <- paste0("/dada_file/","/result_file")
if(!file.exists(wrdir)) dir.create(wrdir, recursive = TRUE)

# specify sample tags from colnames of RPMs
treatment <- c("WT","kro")
reps <- c("r1", "r2", "r3")
tdecay <- c("00", "15", "30", "60", "120", "240", "480") 
#again NO nesting of one label in another hence "00" instead of "0".

mean_RPM <- data.frame("geneID" = rownames(RPMs))
SE_RPM <- data.frame("geneID" = rownames(RPMs))
for(g in treatment){
  for(t in tdecay){
    mean_RPM <- cbind(mean_RPM, rowMeans(RPMs[, cols(df=RPMs, patterns = c(g, t))]))
    names(mean_RPM)[length(mean_RPM)] <- paste0(g, "_", t)
    SE_RPM <- cbind(SE_RPM, apply(X = RPMs[, cols(df=RPMs, patterns = c(g, t))], MARGIN = 1, 
                                  FUN = stats::sd)/sqrt(length(reps)))
    names(SE_RPM)[length(SE_RPM)] <- paste0(g, "_", t)
  }}
mean_RPM[1:20, ]
SE_RPM[1:20, ]
# write output to file
write.table(x = mean_RPM, paste0(wrdir, "/RPM_mean.txt"), sep = "\t")
write.table(x = SE_RPM,   paste0(wrdir, "/RPM_SE.txt"), sep = "\t")

filt1 <- rep(TRUE, 41344) # we will not filter in this example, so the filter value for each gene is set to TRUE.

mT0norm <- data.frame(row.names = rownames(RPMs)[filt1])
for(g in treatment){
  mean_T0reps <- rowMeans(RPMs[filt1, cols(df=RPMs, patterns=c(g, "00"))])
  for(r in reps){
    df <- RPMs[filt1, colnames(RPMs)[cols(df=RPMs, patterns=c(g, r))]]
    df <- df[, 1:length(tdecay)]/mean_T0reps
    mT0norm <- cbind(mT0norm, df)
  }}
write.table(x = mT0norm, file = paste0(wrdir, "/T0 normalized.txt"),  sep = "\t")
tail(mT0norm)

mean_mT0norm <- data.frame(row.names = rownames(mT0norm))
for(g in treatment){
  for(t in tdecay){
    mean_mT0norm <- cbind(mean_mT0norm, rowMeans(mT0norm[, cols(df=mT0norm, patterns=c(g, t))]))
    names(mean_mT0norm)[length(names(mean_mT0norm))] <- paste0(g, "_", t)
  }}

SE_mT0norm <- data.frame(row.names = rownames(mT0norm))
for(g in treatment){
  for(t in tdecay){
    SE_mT0norm <- cbind(SE_mT0norm, apply(X = mT0norm[, cols(df=mT0norm, patterns=c(g, t))], 
                                          MARGIN = 1, FUN = function(x) stats::sd(x)/sqrt(length(reps))))
    names(SE_mT0norm)[length(names(SE_mT0norm))] <- paste0(g, "_", t)
  }}

# write output to file
write.table(x = mean_mT0norm, file = paste0(wrdir, "/T0 normalized_Mean.txt"),  sep = "\t")
write.table(x = SE_mT0norm, file = paste0(wrdir, "/T0 normalized_SE.txt"),  sep = "\t")
tail(mean_mT0norm)
#####################################################################################
##decay factors input
stablegenes <- c( "gene1","gene2","gene3","gene4","gene5",
                  "gene6","gene7","gene8","gene9","gene10",
                  "gene11","gene12","gene13","gene14","gene15",
                  "gene16","gene17","gene18","gene19","gene20",
                  "gene21","gene22","gene23","gene24","gene25",
                  "gene26","gene27","gene28","gene29","gene30")


stabletable <- mean_mT0norm[stablegenes, ]
normFactors <- colMeans(stabletable)
write.table(x <- normFactors, paste0(wrdir, "/Normalziation Decay Factors.txt"), sep = "\t")
normFactors_mean <- matrix(normFactors, nrow = length(tdecay))
normFactors_SE <- matrix(apply(X = stabletable, MARGIN = 2, function(x) stats::sd(x)/sqrt(length(stablegenes))), 
                         nrow = length(tdecay))
t.decay <- c(0,15, 30, 60, 120, 240,480)
rownames(normFactors_mean) <- t.decay
rownames(normFactors_SE) <- t.decay
colnames(normFactors_mean) <- treatment
colnames(normFactors_SE) <- treatment
list(normalizationFactors = normFactors_mean, SE = normFactors_SE)


nF <- vector()
ind <- sapply(names(normFactors), function(x) grep(x, colnames(mT0norm)))
for(i in 1:ncol(ind)){
  nF[ind[, i]] <- colnames(ind)[i]
}
normFactorsM <- t(matrix(rep(normFactors[nF], nrow(mT0norm)), ncol = nrow(mT0norm)))
rm(nF, ind)


mT0norm_2 <- data.frame(mT0norm/normFactorsM, 
                        row.names = rownames(mT0norm))
head(mT0norm_2)
write.table(mT0norm_2, paste0(wrdir, "/T0 normalized and decay factor corrected.txt"), sep = "\t")


mT0norm_2.1 <- reshape2::melt(as.matrix(mT0norm_2), varnames = c("geneID", "variable"))
mT0norm_2.1 <- cbind(mT0norm_2.1, reshape2::colsplit(mT0norm_2.1$variable, "_", names = 
                                                       c("treatment", "t.decay", "rep")))

mT0norm_2.1 <- mT0norm_2.1[, colnames(mT0norm_2.1) !=  "variable"]
mT0norm_2.1 <- mT0norm_2.1[, c(1, 3, 4, 5, 2)]
colnames(mT0norm_2.1) <- c("geneID", "treatment", "t.decay", "rep", "value")
mT0norm_2.1$rep <- gsub("r", "rep", mT0norm_2.1$rep)

mT0norm_2.1$t.decay <- as.numeric(mT0norm_2.1$t.decay)

mT0norm_2.1$treatment <- factor(mT0norm_2.1$treatment,levels = c("WT", "kro"))

mT0norm_2.1$rep <- factor(mT0norm_2.1$rep, levels = paste0("rep",1:3))
write.table(x = mT0norm_2.1,  file = paste0(wrdir, "/ExampleDecayData+stableGenes.txt"), sep = "\t")

tail(mT0norm_2.1)

decay_data  <- mT0norm_2.1 
tail(decay_data)

levels(decay_data$treatment) 

# reorder factor levels if necessary (skipped here):
# decay_data$treatment <- factor(decay_data$treatment, levels = levels(decay_data$treatment)[c(4, 1, 2, 3)]) 
# numbers here refer to the current order position in the factor levels() - the order of the numbers designates the new position

decay_data <- decay_data[order(decay_data$t.decay), ]
decay_data <- decay_data[order(decay_data$rep), ]
decay_data <- decay_data[order(as.numeric(decay_data$treatment)), ] 
decay_data <- decay_data[order(decay_data$geneID), ]
ids <- as.character(unique(decay_data$geneID)) 
tail(decay_data)
write.table(x = decay_data,  file = paste0(wrdir, "/decay_data.txt"), sep = "\t")

nEquivGrp <- if (length(unique(decay_data$treatment)) == 2) {2} else
  if (length(unique(decay_data$treatment)) == 3) {5} else 
    if (length(unique(decay_data$treatment)) == 4) {15}
genoSet <- 1:(length(unique(decay_data$rep)) * length(unique(decay_data$t.decay)))
nTreat <- length(unique(decay_data$treatment))
nSet <- length(genoSet)*nTreat
groups <- groupings(decay_data)
mods <- data.frame(
  "a" = as.vector(sapply(1:nEquivGrp, function(x) {rep(x, nEquivGrp + 1)})),
  "b" = rep(1:(nEquivGrp + 1), nEquivGrp),
  row.names = paste0("mod", 1:(nEquivGrp*(nEquivGrp+1)))
)



group_map(decaydata = decay_data, path = paste0(wrdir, "/Model grouping colormap.pdf"), 
          nEquivGrp = nEquivGrp, groups = groups, mods = mods)


a_bounds <- c(a_low(max(decay_data$t.decay)),
              a_high(min(unique(decay_data$t.decay)[unique(decay_data$t.decay)>0])))

b_bounds <- c(b_low(max(decay_data$t.decay)), 0.075)
a_bounds;b_bounds


###################################################################################
######model test
a <- proc.time()[3]
models <- lapply(X = ids, # alternatively use `ids` here to complete the entire sample data set, but be prepared to wait ~ 10 h. These gene IDs will get passed one at time to the "gene" argument of mod_optimization() and return a list of the results data frame.
                 FUN = mod_optimization,
                 data = decay_data, 
                 group = groups, 
                 mod = mods,
                 alpha_bounds = a_bounds, 
                 beta_bounds = b_bounds,
                 models = rownames(mods), 
                 path = paste0(wrdir, "/modeling_results"),
                 file_only = FALSE) 
names(models) <- ids
b <- proc.time()[3]
(b-a)/60/length(ids) # gives you average min per gene



models <- lapply( paste0( wrdir, "/modeling_results/", ids, "_results.txt"), read.delim, header = TRUE )   
names(models) <- ids

##########################################################################################
##model selection
## -----------------------------------------------------------------------------
results <- t(sapply(models, function(x) x[x[, "AICc"] == min(x[, "AICc"]), ]))
results <- as.data.frame(results)

results[, 1:2] <- sapply(as.data.frame(results[, 1:2]), function(x) as.character(unlist(x)))
results[, -c(1,2)] <- sapply(results[, -c(1,2)], unlist)
write.table(results, file = paste0(wrdir,"/best model results.txt"), sep = "\t")
results <- read.delim(paste0(wrdir,"/best model results.txt"))

##
results <- read.delim(paste0(wrdir,"/best model results.txt"))
results$alpha_grp <- mods[as.character(results$mod), "a"]
results$beta_grp <- mods[as.character(results$mod), "b"]
results$mod <- as.numeric(gsub("mod", "", as.character(results$mod)))

results$alphaPattern <- sapply(rownames(results), function(x) {
  paste0(gsub("alpha_", "", colnames(results)[3:(2+nTreat)][order(round(results[x, 3:(2+nTreat)], 4))]), collapse = "<=")
})
results$alphaPattern <- paste0(results$alpha_grp, "_", results$alphaPattern)
results$betaPattern <- sapply(rownames(results), function(x){
  paste0(gsub("beta_", "", colnames(results)[(3+nTreat):(2+2*nTreat)][order(round(results[x, (3+nTreat):(2+2*nTreat)], 4))]), collapse = "<=")
})
results$betaPattern <- paste0(results$beta_grp, "_", results$betaPattern)

results <- results[order(rownames(results)), ]
results <- results[order(results$beta_grp), ]
results <- results[order(results$alphaPattern), ]
results <- results[order(results$alpha_grp), ]

results$alphaPattern <- factor(results$alphaPattern, levels = as.character(unique(results$alphaPattern)))

results <- data.frame(results[, 3:(2*nTreat+3), 2], results[, c("AICc", "alpha_grp", "beta_grp", "alphaPattern", "betaPattern")])
results$nEqMods <- sapply(min_mods[rownames(results)], length)
results$nEqAgp <- sapply(min_alpha_mods[rownames(results)], length)

write.table(results, paste0(wrdir,"/alphas+betas+mods+grps+patterns+relABs.txt"), sep = "\t")