library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library("UpSetR")
library(corrplot)

#reading output from derfinder
table3 <- read.table("GBM_Output2.txt")


#adding rownames as chr:start-end
colnames(table3) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(table3))
rownames(table3) <- paste(table3$contig,":",table3$start,"-",table3$end,sep="")

#Just to remove the decimals in the count from derfinder
count3 <- round(as.data.frame(table3[,4:27]))

#remove some sufixes from colnames
colnames(count3) <- sub(".RNA.NovaSeq_.*","",colnames(count3))

#Create object for edgeR
y4 <- DGEList(counts=count3,genes=row.names(count3))

#Keep rows with more than 1 read in more than 15 samples
keep <- rowSums(cpm(y4)>1) >= 15
y45 <- y4[keep,]

#normalization with TMM
y45 <- calcNormFactors(y45, method="TMM")
y45 <- estimateCommonDisp(y45)
extract cpm values
cpms45 <- cpm(y45)

#plot variability
dist_mat45 <- dist(t(cpms45), method = 'euclidean')
hclust_avg45 <- hclust(dist_mat45, method = 'complete')
plot(hclust_avg45)
#Two ways of doing MDS plot
plotMDS(y4[,-19])
cmds45 <- cmdscale(dist_mat45, eig = TRUE, k = 2)
x45 <- cmds45$points[, 1]
y45 <- cmds45$points[, 2]
ggplot() + geom_point(data = as.data.frame(cmds45$eig) , mapping = aes(x = x45, y = y45), color = c(rep("red",3),rep("green",4),rep("blue",5),rep("black",2),rep("pink",3),rep("brown",2),"orange",rep("yellow",2),"white","salmon"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(count3) ,cmds45$eig), mapping = aes(x=x45, y=y45 , label = colnames(count3))  ) + labs(title = "figure5:  MDS configuration of Japan Prefectures with labels")


#organizing by clusters 
count_CDEFG2 <- count_CDEFG2[,c(2,4,6,17,9,10,11,16,3,5,7,12:15,8,18:23,1)] 
designG3 <- data.frame(row.names= colnames(count_CDEFG2),sample= colnames(count_CDEFG2),condition=c(rep("Progression",8),rep("NoProgression",15)))
designG3
designG3$condition <- factor(designG3$condition,levels=c("NoProgression","Progression"))
design.matrixG3 <- model.matrix(~ condition,data=designG3)
#filtering by 10 counts at least in 4 samples in NoProgression and 6 samples in Progression
keepCDEFG2 <- which((rowSums(count_CDEFG2[,1:8] > 10) >= 4 ) | (rowSums(count_CDEFG2[,9:23] > 10) >= 6))
count_CDEFG_filt2 <- count_CDEFG2[names(keepCDEFG2),]
yCDEFG2 <- DGEList(counts=count_CDEFG_filt2,genes=row.names(count_CDEFG_filt2))
yCDEFG2 <- calcNormFactors(yCDEFG2, method="TMM")
yCDEFG2 <- estimateCommonDisp(yCDEFG2)
plotMDS(yCDEFG2)
#differential expression
yCDEFG2<- estimateGLMTrendedDisp(yCDEFG2,design.matrixG3)
yCDEFG2 <- estimateGLMTagwiseDisp(yCDEFG2,design.matrixG3)
fitcG3 <- glmFit(yCDEFG2,design.matrixG3)
rescGG3 <- glmLRT(fitcG3,coef=2)
resultsG3 <- topTags((rescGG3),n=nrow(yCDEFG2$rows))
topTags(rescGG3)
subset(resultsG$table,resultsG$table$FDR < 0.05)

#print out results
write.table (resultstG3, file="results_GBM_groupsALL.tsv",col.names = T,row.names=T,sep="\t")

#just convert sample names to original
convert_list2$original_name =c ("BS20-105C","BS20-113D","BS20-105E","BS20-110F","BS20-105F","BS20-113F","BS20-111F","BS20-110E","BS20-113C","BS20-111D","BS20-113E","BS20-110C","BS20-117C","BS20-110D","BS20-123D","BS20-111E","BS20-114C","BS21-126C","BS20-112D","BS21-126D","BS21-128D","BS20-112F","BS20-113G")
convert_list2
colnames(FDR_sigG3)[24:46] <- convert_list2$original_name
cpms_FDR3 <- FDR_sigG3[,24:46]
head(cpms_FDR3)
#heatmap
heatmap1 <- heatmap.2(log10(as.matrix(cpms_FDR3+0.5/rowMeans(cpms_FDR3+0.5))),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(149),ColSideColors = c(rep("red",8),rep("green",15)),margins=c(8,8), breaks=c(seq(-7,-1,length=50),seq(-0.9,0.9,length=50), seq(1,7,length=50)),main="Total FDR < 0.05",labRow=F)
sessionInfo()

#ranking analysis 
#Getting the cpms of all significant DE regions
trans2 <- as.data.frame(t(cpms_FDR3))
#adding column with the conditions 1 for no progression and 0 for progression
trans2$diag <- c(rep(1,8),rep(0,15))
#using correlation bigger than 75% between the condition and each DE region
trans_red2 <- trans2[,names(cor(trans2)["diag",which(cor(trans2)["diag",] > 0.75)])]
