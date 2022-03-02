tablealz <- read.table("alzheimer.Output.txt")
spikeina <- read.table("spikein_counts.tzt",header=T)
head(table2)

colnames(tablealz) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(tablealz))
rownames(tablealz) <- paste(tablealz$contig,":",tablealz$start,"-",tablealz$end,sep="")
head(tablealz)
counta <- round(as.data.frame(tablealz[,4:23]))
colnames(spikein)
head(spikein[,colnames(count)])
head(count)
count2a <- rbind(counta,spikeina[,colnames(counta)])
head(count2a)
head(gsub(".RNA.NovaSeq_S\d+","",colnames(count2a)))
colnames(count2a) <- sub(".RNA.NovaSeq_.*","",colnames(count2a))
rownames(count2a) = c(rownames(tablealz), spikeina$spikein)
ya <- DGEList(counts=count2a,genes=row.names(count2a))
keep <- rowSums(cpm(ya)>1) >= 12
y2a <- ya[keep,]

y2a <- calcNormFactors(y2a, method="TMM")
y2a <- estimateCommonDisp(y2a)
rownames(y2a)

spike_normA <- cpm(y2a)[922:960,]
par(cex=0.5)
corrplot(cor(spike_normA),method= "number")
par(cex=1)
plotMDS(y2a[1:921,])
cpmsA <- cpm(y2a)[1:921,]
dist_matA <- dist(t(cpmsA), method = 'euclidean')
hclust_avgA <- hclust(dist_matA, method = 'average')
dend <- as.dendrogram(hclust_avg)
dend <- color_labels(dend,k=3)
dend
length(colnames(cpmsA))
codes <- c("red","green","red","green",rep("red",2),"green",rep("red",3),"green","red","green","red","green","red","green",rep("red",2))
length(codes)
label_dolors(dend) <- 
str(hclust_avgA)
str(dend)
tablealz2 <- read.table("alzheimer2.Output.txt")
colnames(tablealz2) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(tablealz2))
rownames(tablealz2) <- paste(tablealz2$contig,":",tablealz2$start,"-",tablealz2$end,sep="")
colnames(tablealz2)
counta2 <- round(as.data.frame(tablealz2[,4:22]))
#counta2 <- 
colnames(counta2) <- sub(".RNA.NovaSeq_.*","",colnames(counta2))
counta2<-counta2[,c(2,4,7,11,13,15,17,1,3,5,6,8,9,10,12,14,16,18,19)] 
colnames(counta2)
groupsA <- c(rep(1,7),rep(2,12))
batchesA <- c(1,2,2,3,3,1,1,1,1,2,2,2,2,3,3,2,2,3,1)
correctA <- ComBat_seq(counts = as.matrix(counta2),batch=batchesA,group=groupsA)
length(batchesA)
head(counta2)

ya2 <- DGEList(counts=counta2,genes=row.names(counta2))
keep <- rowSums(cpm(ya2)>1) >= 12
ya3 <- ya2[keep,]
#ya3 <- ya19)]
ya3 <- calcNormFactors(ya3, method="TMM")
ya3 <- estimateCommonDisp(ya3)
cpmsA2 <- cpm(ya3)
dist_matA2<- dist(t(cpmsA2), method = 'euclidean')
hclust_avgA2 <- hclust(dist_matA2, method = 'average')
plot(hclust_avgA2)
colnames(ya3)
codes <- c("green","green","green","green",rep("green",2),"green",rep("red",3),"red","red","red","red","red","red","red",rep("red",2))
dend <- as.dendrogram(hclust_avgA2)
labels_colors(dend) <- codes[order.dendrogram(dend)]
pdf("Clustering_data.pdf")

plot(dend,main="All Samples Clustering - All features")
mdsA <- plotMDS(ya3)
barplot(ya3$samples$lib.size)
designA
data.xA <- mdsA$x
data.yA <- mdsA$y
data.z <- mdsef$y
pdf("mds_DEF_multiple.pdf")
data.mergedA <- data.frame(data.xA,data.yA)
ggplot() + geom_point(data = as.data.frame(data.mergedA) , mapping = aes(x = data.xA, y = data.yA), color = c("green","green","green","green",rep("green",2),"green",rep("red",3),"red","red","red","red","red","red","red",rep("red",2)), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpmsA2) ,data.mergedA), mapping = aes(x=data.xA, y=data.yA , label = colnames(cpmsA2))  ) + labs(title = "MDS Plot - All Samples - All features",x="component 1",y="component 2")

tablealz3 <- read.table("alzheimer3.Output.txt")
colnames(tablealz3) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(tablealz3))
rownames(tablealz3) <- paste(tablealz3$contig,":",tablealz3$start,"-",tablealz3$end,sep="")
colnames(tablealz3)
counta3 <- round(as.data.frame(tablealz3[,4:18]))
colnames(counta3) <- sub(".RNA.NovaSeq_.*","",colnames(counta3))
colnames(counta3)
counta3 <- counta3[,-15]
ya5 <- DGEList(counts=counta3,genes=row.names(counta3))
keepAs <- rowSums(cpm(ya5)>1) >= 6
ya6 <- ya5[keepAs,]
plotMDS(ya6)
mdsA2 <- plotMDS(ya6)
data.xA2 <- mdsA2$x
data.yA2 <- mdsA2$y
data.mergedA2 <- data.frame(data.xA2,data.yA2)
data.mergedA2 
colnames(cpm(ya6))
ggplot() + geom_point(data = as.data.frame(data.mergedA2) , mapping = aes(x = data.xA2, y = data.yA2), color = c(rep("green",5),rep("red",9)), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(ya6)) ,data.mergedA2), mapping = aes(x=data.xA2, y=data.yA2 , label = colnames(cpm(ya6)))  ) + labs(title = "MDS Plot - No BNM7,8,9,10 - All features",x="component 1",y="component 2")
color_sex1 <- c(rep("gold",2),rep("darkorchid",6),rep("gold",2),"darkorchid","gold",rep("darkorchid",2))
ggplot() + geom_point(data = as.data.frame(data.mergedA2) , mapping = aes(x = data.xA2, y = data.yA2), color = c(rep("green",5),rep("red",9)), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(ya6)) ,data.mergedA2), mapping = aes(x=data.xA2, y=data.yA2 , label = colnames(cpm(ya6)))  ) + labs(title = "MDS Plot - No BNM7,8,9,10 - All features",x="component 1",y="component 2")
ggplot() + geom_point(data = as.data.frame(data.mergedA2) , mapping = aes(x = data.xA2, y = data.yA2), color = color_sex1, alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(ya6)) ,data.mergedA2), mapping = aes(x=data.xA2, y=data.yA2 , label = colnames(cpm(ya6)))  ) + labs(title = "MDS Plot - No BNM7,8,9,10,18 - All features",x="component 1",y="component 2")
dist_matA3<- dist(t(cpm(ya6)), method = 'euclidean')
hclust_avgA3 <- hclust(dist_matA3, method = 'average')
plot(hclust_avgA3)
dend4 <- as.dendrogram(hclust_avgA3)
colnames(cpm(ya6))
codes3<- c(rep("green",5),rep("red",9))
pdf("sex_color_coding.pdf")
labels_colors(dend4) <- code3[order.dendrogram(dend4)]
plot(dend4,main="No BNM7,8,9,10,18 Clustering - All features")
labels_colors(dend4) <- color_sex1[order.dendrogram(dend4)]
plot(dend4,main="No BNM7,8,9,10,18 Clustering - All features")
colnames(counta3)
designa <- data.frame(row.names= colnames(counta3),sample= colnames(counta3),condition=c(rep("AD",5),rep("Control",9)))
designa$condition <- factor(designa$condition,levels=c("AD","Control"))
design.matrixa <- model.matrix(~ condition,data=designa)
keepRa <- which((rowSums(ya5$counts[,1:5] > 10) >=3) | (rowSums(ya5$counts[,6:14] > 10) >=5))
ya7 <- ya5[names(keepRa),]
ya7 <- calcNormFactors(ya7,method="TMM")
ya7 <- estimateGLMCommonDisp(ya7,design.matrixa)
barplot(ya6$samples$lib.size)
ya7 <- estimateGLMTrendedDisp(ya7,design.matrixa)
ya7 <- estimateGLMTagwiseDisp(ya7,design.matrixa)
fita <- glmFit(ya7,design.matrixa)
resa <- glmLRT(fita,coef=2)
topTags((resa))
totala <- topTags(resa,n = nrow(ya6$counts))
countsat <- cbind(round(ya7$counts),cpm(ya7))
resultsa <- cbind(countsat,totala$table[rownames(countsat),])
resultsla2 <- cbind(resultsa,tablealz3[rownames(resultsa),19:24])
head(resultsla2[order(resultsla2$PValue),],20)
dim(resultsla2)
fdra <- subset(resultsla2[,15:28],resultsla2$FDR < 0.05)
fdrT <- subset(resultsla2,resultsla2$FDR < 0.05)
write.table(fdrT,file="fdr_results.csv",sep="\t",quote = F)
fdr_b <- read.table("fdr_binary.csv",header=T)
pdf("upset_plots.pdf")
upset(fdr_b,nset=8,keep.order=T,nintersects = 62,set_size.show=F)
fdrT
pvalueT <- subset(resultsla2,resultsla2$PValue < 0.05)
write.table(pvalueT,file="pvalue_results.csv",sep="\t",quote = F)
pvalue_b <- read.table("pvalue_binary.csv",header=T)
upset(pvalue_b,nset=8,keep.order=T,nintersects = 62,set_size.show=F)
dev.off()
colnames(fdrT)
dim(fdra)
heatmap.2(as.matrix(log2((fdra+1)/rowMeans(fdra+1))),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(299),ColSideColors = c(rep("green",5),rep("red",9)),margins=c(8,8), breaks=c(seq(-5,-1,length=100),seq(-0.9,0.9,length=100), seq(1,5,length=100)),main="Heatmap FDR < 0.05")
write.table(resultsla2,"alzheimer_comparision_derfinder.tsv",sep="\t")
groupsA <- c(rep(1,5),rep(2,9))
batchesA <- c(1,2,2,1,1,1,1,rep(2,6),1)
length(batchesA)
correctA <- ComBat_seq(counts = as.matrix(counta3),batch=batchesA,group=groupsA)
yaA <- DGEList(counts=correctA,genes=row.names(correctA))
keepAsC <- rowSums(cpm(yaA)>1) >= 6
yaA2 <- yaA[keepAsC,]
mdsAC <- plotMDS(yaA2)
data.xAC <- mdsAC$x
data.yAC <- mdsAC$y
data.mergedAC <- data.frame(data.xAC,data.yAC)
ggplot() + geom_point(data = as.data.frame(data.mergedAC) , mapping = aes(x = data.xAC, y = data.yAC), color = c(rep("green",5),rep("red",9)), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(yaA2)) ,data.mergedAC), mapping = aes(x=data.xAC, y=data.yAC , label = colnames(cpm(yaA2)))  ) + labs(title = "MDS Plot - No BNM7,8,9,10 - All features",x="component 1",y="component 2")




excerptmirna <- read.table("exceRpt_miRNA_ReadCounts.txt")
head(excerptmirna)
colnames(excerptmirna) <- sub(".RNA.NovaSeq_.*","",colnames(excerptmirna))
head(excerptmirna)
excerptmirna2 <- excerptmirna[,-19]
excerptmirna2 <- excerptmirna2[,c(2,4,7,11,13,15,17,1,3,5,6,8,9,10,12,14,16,18,19)]
exmirna <- DGEList(counts=excerptmirna2,genes=row.names(excerptmirna2))
keepe <- rowSums(cpm(exmirna)>1) >= 12
exmirna2 <- exmirna[keepe,]
dist_exmirna <- dist(t(cpm(exmirna2)), method = 'euclidean')
hclust_exmirna <- hclust(dist_exmirna, method = 'complete')
plot(hclust_exmirna)
dend2 <- as.dendrogram(hclust_exmirna)
codes2<- c(rep("green",7),rep("red",12))
labels_colors(dend2) <- codes2[order.dendrogram(dend2)]
plot(dend2,main="All Samples clustering-Only miRNAs")
mdsExmir <- plotMDS(exmirna2)
data.xEx <- mdsExmir$x
data.yEx <- mdsExmir$y
data.mergedExmir <- data.frame(data.xEx,data.yEx)
ggplot() + geom_point(data = as.data.frame(data.mergedExmir) , mapping = aes(x = data.xEx , y = data.yEx), color = codes2, alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(exmirna2)) ,data.mergedExmir), mapping = aes(x=data.xEx, y=data.yEx , label = colnames(cpm(exmirna2)))  ) + labs(title = "MDS all samples-miRNA",x="Component 1", y= "Component 2")
head(excerptmirna3)
excerptmirna3 <- excerptmirna2[,-c(4,5,14,15)]
excerptmirna3 <- excerptmirna3[,-15]
exmirna4 <- DGEList(counts=excerptmirna3,genes=row.names(excerptmirna3))
keepe2 <- rowSums(cpm(exmirna4)>1) >= 7
exmirna5 <- exmirna4[keepe2,]
dist_exmirna2 <- dist(t(cpm(exmirna5)), method = 'euclidean')
hclust_exmirna2 <- hclust(dist_exmirna2, method = 'complete')
plot(hclust_exmirna2)
dend3 <- as.dendrogram(hclust_exmirna2)
codes3<- c(rep("green",5),rep("red",10))
labels_colors(dend3) <- codes3[order.dendrogram(dend3)]
plot(dend3,main="No Samples BNM7,8,9,10,18 clustering-Only miRNAs")
labels_colors(dend3) <- color_sex1[order.dendrogram(dend3)]
plot(dend3,main="No Samples BNM7,8,9,10,18 clustering-Only miRNAs")
plotMDS(exmirna2)
dev.off()
plot(dend3,main="No Samples BNM7,8,9,10,18 clustering-Only miRNAs")
mdsExmir2 <- plotMDS(exmirna5)
mdsExmir3 <- plotMDS(exmirna5,dim.plot = c(1,3))
data.xEx2 <- mdsExmir2$x
data.yEx2 <- mdsExmir2$y
data.zEx2 <- mdsExmir3$y
data.mergedExmir2 <- data.frame(data.xEx2,data.yEx2,data.zEx2)

ggplot() + geom_point(data = as.data.frame(data.mergedExmir2) , mapping = aes(x = data.xEx2 , y = data.yEx2), color = color_sex1, alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(exmirna5)) ,data.mergedExmir2), mapping = aes(x=data.xEx2, y=data.yEx2 , label = colnames(cpm(exmirna5)))  ) + labs(title =  "MDS Plot - No BNM7,8,9,10,18 - All features",x="component 1",y="component 2")
ggplot() + geom_point(data = as.data.frame(data.mergedExmir2) , mapping = aes(x = data.xEx2 , y = data.zEx2), color = codes3, alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm(exmirna5)) ,data.mergedExmir2), mapping = aes(x=data.xEx2, y=data.zEx2 , label = colnames(cpm(exmirna5)))  ) + labs(title = "MDMS Alzheimer 1x3")

head(exmirna4)
keepFa <- which((rowSums(exmirna4$counts[,1:5] > 10) >=3) | (rowSums(exmirna4$counts[,6:14] > 10) >=5))
head(keepFa)
dim(exmirna5$counts)
exmirna6 <- exmirna4[names(keepFa),]
dim(exmirna6$counts)
exmirna6 <-calcNormFactors(exmirna6,method="TMM")
exmirna6<- estimateGLMCommonDisp(exmirna6,design.matrixa)
barplot(ya6$samples$lib.size)
exmirna6 <- estimateGLMTrendedDisp(exmirna6,design.matrixa)
exmirna6 <- estimateGLMTagwiseDisp(exmirna6,design.matrixa)
fitea <- glmFit(exmirna6,design.matrixa)
resea <- glmLRT(fitea,coef=2)
topTags((resea))
totalea <- topTags(resea,n = nrow(exmirna6$counts))
countseat <- cbind(round(exmirna6$counts),cpm(exmirna6))
resultsea <- cbind(countseat,totalea$table[rownames(countseat),])
#resultsla2 <- cbind(resultsa,tablealz3[rownames(resultsa),19:24])
head(resultsea[order(resultsea$FDR),],20)
head(resultsea)
resultsea["hsa-miR-203b",]
write.table(resultsea,"alzheimer_comparision_mirnaonly.tsv",sep="\t")
dev.off()

