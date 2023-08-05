#----------------------------------------------------------
# TODO 统计每个基因在两种方法之间的兰德指数
#----------------------------------------------------------
library(flexclust)
library(dplyr)
library(tidyr)

stageArray=data.frame(
    stage=c('test'),
    meanVal=c(2),
    minVal=c(3),
    maxVal=c(4)
)
for (stage in c("0DPA", "4DPA", "8DPA", "12DPA", "16DPA", "20DPA"))
{
    RIarray <- c()
    FPKM_BiasType <- read.table(paste(stage, "/", stage, "_BiasType_Sample_FPKM.txt", sep = ""), header = T, row.names = 1, sep = "\t")
    reads_BiasType <- read.table(paste(stage, "/", stage, "_BiasType_Sample_readCount.txt", sep = ""), header = T, row.names = 1, sep = "\t")
    commonSample <- intersect(colnames(reads_BiasType), colnames(FPKM_BiasType))
    FPKM_BiasType <- FPKM_BiasType %>% select(as.vector(commonSample))
    reads_BiasType <- reads_BiasType %>% select(as.vector(commonSample))
    for (genePair in rownames(reads_BiasType)) {
        readCount <- reads_BiasType[genePair, ]
        FpkmCount <- FPKM_BiasType[genePair, ]
        RIindex <- randIndex(table(readCount, FpkmCount), correct = F)
        RIarray <- c(RIarray, RIindex)
    }
    print(paste(stage, mean(RIarray), min(RIarray), max(RIarray), sep = "--"))
    stageArray=rbind(stageArray,list(stage, mean(RIarray), min(RIarray), max(RIarray)))
}
print(stageArray)
write.table(stageArray,"geneRandIndex.txt",col.names =T,row.names = F,quote = F,sep="\t")