library(edgeR)
tableall <- read.table("Output.txt")
colnames(tableall) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(tableall))
rownames(tableall) <- paste(tableall$contig,":",tableall$start,"-",tableall$end,sep="")
countl <- round(as.data.frame(tableall[,4:13]))
count100 <- countl[,c(1,2)]
count100red <- count100[names(which((rowSums(count100 > 10) >=1))),] 
plot(log2(count100red))
summary(lm(count100red))
y100 <- DGEList(counts=count100red,genes=row.names(count100red))
y100 <- calcNormFactors(y100, method="TMM")
y100n <- cpm(y100)
plot(log2(y100n))


head(tableall)
count_r <- tableall[,c(4,6,7,8,9,10,12,13)]
head(count_r)
red_final_der <- count_r
head(red_final_der)



red_final_der <- red_final_der[,c(2,3,8,5,6,7,1,4)]
keepNonxRd <- which((rowSums(red_final_der[,c(1:3)] > 10) >= 2 ) | (rowSums(red_final_der[,c(4:6)] > 10) >= 2))
head(red_final_mirna_filt)
red_final_der_filt <- red_final_der[names(keepNonxRd),]
length(keepNonxR)
y_der <- DGEList(counts=red_final_der_filt,genes=row.names(red_final_der_filt))
y_der <- estimateCommonDisp(y_der,method="TMM")
plotMDS(y_der)
y_der2 <- y_der[,-c(7,8)]
plotMDS(y_der2)
design2
y_der2<- estimateGLMTrendedDisp(y_der2,design.matrix)
y_der2 <- estimateGLMTagwiseDisp(y_der2,design.matrix)
fit_der <- glmFit(y_der2,design.matrix)
resder <- glmLRT(fit_der,coef=2)
head(resmirna)
results_der <- topTags((resder),n=nrow(y_der2$rows))
summary(results_der$table$FDR)
rownamesPvalued <- rownames(subset(results_der$table,results_der$table$PValue < 0.05))
head(rownamesPvalue)
y_der3 <- y_der[rownamesPvalued,]
plotMDS(y_der3)
countsder <- cbind(y_der2$counts,cpm(y_der2))
head(countsl)
head(totall)
resultsder <- cbind(countsder,results_der$table[rownames(countsder),])

resultsderT <- cbind(resultsder,tableall[rownames(resultsder),14:19])
write.table (resultsdercT, file="derfinder_Notcorrected.tsv",col.names = T,row.names=T,sep="\t")
###correction 
library(sva)
head(y_der$counts)
groupsl <- c(rep(1,3),rep(2,5))
batchesl <- c(1,2,2,2,1,1,1,2)

correctl <- ComBat_seq(counts = as.matrix(y_der$counts),batch=batchesl,group=groupsl)
head(correctl)
y_derc <- DGEList(counts=correctl,genes=row.names(correctl))
plotMDS(y_derc)
design2 <- data.frame(row.names= colnames(correctl),sample= colnames(correctl),condition=c(rep("Responder",3),rep("NonResponder",5)))
design2$condition <- factor(design2$condition,levels=c("Responder","NonResponder"))
design.matrix2 <- model.matrix(~ condition,data=design2)
y_derc<- estimateGLMTrendedDisp(y_derc,design.matrix2)
y_derc <- estimateGLMTagwiseDisp(y_derc,design.matrix2)
fit_derc <- glmFit(y_derc,design.matrix2)
resderc <- glmLRT(fit_derc,coef=2)
head(resmirna)
results_derc <- topTags((resderc),n=nrow(y_derc$rows))
dim(subset(results_derc$table, results_derc$table$FDR < 0.05))
rownamesFDRc <- rownames(subset(results_derc$table,results_derc$table$FDR < 0.05))
head(rownamesPvalue)
y_derc3 <- y_derc[rownamesFDRc,]
pmds4 <- plotMDS(y_derc3)
cpmderc <- cpm(y_derc)
head(pmds)
xcl4 <- pmds4$eigen.vectors[,1]
ycl4 <- pmds4$eigen.vectors[, 2]
ggplot() + geom_point(data = as.data.frame(pmds4$eigen.values) , mapping = aes(x = xcl4, y = ycl4), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpmderc) ,pmds4$eigen.values), mapping = aes(x=xcl4, y=ycl4 , label = colnames(cpmderc)  )) + labs(title = "MDS_definder_corrected_lung")
countsderc <- cbind(y_derc$counts,cpm(y_derc))
head(countsl)
head(totall)
resultsderc <- cbind(countsderc,results_derc$table[rownames(countsderc),])
head(resultsmrnac)
library(EnhancedVolcano)
colnames(resultsderc)
resultsdercV <- resultsderc[,17:22]
resultsdercV$diffexpressed <- "NO"
resultsdercV$diffexpressed[resultsdercV$logFC > 2 & resultsdercV$FDR < 0.05] <- "UP"
resultsdercV$diffexpressed[resultsdercV$logFC < -2 & resultsdercV$FDR < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
ggplot(data=resultsdercV, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) 
fdrderc <- subset(resultsderc[,9:16],resultsderc$FDR < 0.01)
dim(fdrderc)
colnames(fdrmrna)
heatmap.2(as.matrix(log2((fdrderc+1)/rowMeans(fdrderc+1))),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),labRow=F,main="Derfinder corrected DE")
colnames(tableall)
resultsdercT <- cbind(resultsderc,tableall[rownames(resultsderc),14:19])
write.table (resultsdercT, file="derfinder_corrected.tsv",col.names = T,row.names=T,sep="\t")


#########################
mirna_list <- list.files(path=".",recursive = T, pattern="miRNAmature_sense.txt")
mirna_list
final <- read.table("L100_100ng.fastq/readCounts_miRNAmature_sense.txt",header=T)
final <- final[,-c(2,3,5)]
head(temp)
for (i in 2:length(mirna_list)) {
  temp <- read.table(mirna_list[i],header=T)
  temp <- temp[,-c(2,3,5)]
  final <- merge(final,temp,by="ReferenceID",all=T)
  
  
}
head(final)
colnames(final)[2:11] <- colnames(tableall)[4:13]
head(final)
final[is.na(final)] <- 0
red_final_mirna <- final[,-c(3,9)]
head(red_final_mirna)
rownames(red_final_mirna) <- red_final_mirna$ReferenceID
red_final_mirna$ReferenceID <- NULL
head(red_final_mirna)
red_final_mirna <- red_final_mirna[,c(2,3,8,5,6,7,1,4)]
red_fnal_mirna_filt <- red_final_mirna[names(which(rowMeans(red_final_mirna) > 10)),]
design <- data.frame(row.names= colnames(red_final_mirna),sample= colnames(red_final_mirna),condition=c(rep("Responder",3),rep("NonResponder",3)))
design$condition <- factor(design$condition,levels=c("NonResponder","Responder"))
design.matrix <- model.matrix(~ condition,data=design)
keepNonxR <- which((rowSums(red_final_mirna[,c(1:3)] > 10) >= 2 ) | (rowSums(red_final_mirna[,c(4:6)] > 10) >= 2))
head(red_final_mirna_filt)
red_final_mirna_filt <- red_final_mirna[names(keepNonxR),]
length(keepNonxR)
y_mirna <- DGEList(counts=red_final_mirna_filt,genes=row.names(red_final_mirna_filt))
y_mirna <- estimateCommonDisp(y_mirna,method="TMM")
plotMDS(y_mirna)
y_mirna2 <- y_mirna[,-c(7,8)]
y_mirna2<- estimateGLMTrendedDisp(y_mirna2,design.matrix)
y_mirna2 <- estimateGLMTagwiseDisp(y_mirna2,design.matrix)
fit_mirna <- glmFit(y_mirna2,design.matrix)
resmirna <- glmLRT(fit_mirna,coef=2)
head(resmirna)
results_mirna <- topTags((resmirna),n=nrow(y_mirna2$rows))
summary(results_mirna$table$PValue)
rownamesPvalue <- rownames(subset(results_mirna$table,results_mirna$table$PValue < 0.05))
head(rownamesPvalue)
y_mirna3 <- y_mirna[rownamesPvalue,]
plotMDS(y_mirna3)

cpmymirna <- cpm(y_mirna)
heatmap.2(log2(cpmymirna+1/(rowMeans(cpmymirna) +1)))
dim(cpmyn)
countsmirna <- cbind(y_mirna2$counts,cpm(y_mirna2))
head(countsl)
head(totall)
resultsmirna <- cbind(countsmirna,results_mirna$table[rownames(countsmirna),])
write.table (resultsmirna, file="excerpt_mirna_NOTcorrected.tsv",col.names = T,row.names=T,sep="\t")
library(gplots)
pdf(heatmap.2(cpmyn),file=heatmap.pdf)
#####correction
groupsl <- c(rep(1,3),rep(2,5))
batchesl <- c(1,2,2,2,1,1,1,2)
correctm <- ComBat_seq(counts = as.matrix(y_mirna$counts),batch=batchesl,group=groupsl)
head(correctl)
y_mirnac <- DGEList(counts=correctm,genes=row.names(correctm))
plotMDS(y_mirnac)
#design2mi <- data.frame(row.names= colnames(correctm),sample= colnames(correctm),condition=c(rep("Responder",3),rep("NonResponder",5)))
#design2$condition <- factor(design2$condition,levels=c("NonResponder","Responder"))
#design.matrix2 <- model.matrix(~ condition,data=design2)
y_mirnac<- estimateGLMTrendedDisp(y_mirnac,design.matrix2)
y_mirnac <- estimateGLMTagwiseDisp(y_mirnac,design.matrix2)
fit_mirnac <- glmFit(y_mirnac,design.matrix2)
resmirnac <- glmLRT(fit_mirnac,coef=2)
head(resmirna)

results_mirnac <- topTags((resmirnac),n=nrow(y_mirnac$rows))
dim(subset(results_mirnac$table, results_mirnac$table$FDR < 0.05))
rownamesFDRmirnac <- rownames(subset(results_mirnac$table,results_mirnac$table$FDR < 0.05))
head(rownamesPvalue)
y_mirnac3 <- y_mirnac[rownamesFDRmirnac,]
plotMDS(y_mirnac3)
pmds3 <- plotMDS(y_mirnac3)
head(pmds)
xcl3 <- pmds3$eigen.vectors[,1]
ycl3 <- pmds3$eigen.vectors[, 2]
cpm_mirnac <- cpm(y_mirnac)
ggplot() + geom_point(data = as.data.frame(pmds3$eigen.values) , mapping = aes(x = xcl3, y = ycl3), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mirnac) ,pmds3$eigen.values), mapping = aes(x=xcl3, y=ycl3 , label = colnames(cpm_mirnac)  )) + labs(title = "MDS_mirna_lung")
countsmirna <- cbind(y_mirnac$counts,cpm(y_mirnac))
head(countsl)
head(totall)
resultsmirnac <- cbind(countsmirna,results_mirnac$table[rownames(countsmirna),])
head(resultsmrnac)
library(EnhancedVolcano)

resultsmirnacV <- resultsmirnac[,17:22]
resultsmirnacV$diffexpressed <- "NO"
resultsmirnacV$diffexpressed[resultsmirnacV$logFC > 2 & resultsmirnacV$FDR < 0.05] <- "UP"
resultsmirnacV$diffexpressed[resultsmirnacV$logFC < -2 & resultsmirnacV$FDR < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
ggplot(data=resultsmirnacV, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) 
fdrmirna <- subset(resultsmirnac[,9:16],resultsmirnac$FDR < 0.05)
colnames(fdrmrna)
rownames(fdrmirna) <- gsub(":MIMAT\\d+:Homo:sapiens:miR.+","",rownames(fdrmirna))
heatmap.2(as.matrix(log2((fdrmirna+1)/rowMeans(fdrmirna+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,12), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="miRNA Excerpt corrected DE")
write.table (resultsmirnac, file="excerpt_mirna_corrected.tsv",col.names = T,row.names=T,sep="\t")
#######mRNA
mrna_list <- list.files(path=".",recursive = T, pattern="gencode_sense_geneLevel.txt")
mrna_list
finalm <- read.table("L100_100ng.fastq/readCounts_gencode_sense_geneLevel.txt",header=T)
finalm <- finalm[,-c(2,3,5)]
head(finalm)
head(temp)
for (i in 2:length(mrna_list)) {
  temp <- read.table(mrna_list[i],header=T)
  temp <- temp[,-c(2,3,5)]
  finalm <- merge(finalm,temp,by="GeneSymbol",all=T)
  
  
}

colnames(finalm)[2:11] <- colnames(tableall)[4:13]
head(finalm)
finalm[is.na(finalm)] <- 0
red_final_mrna <- finalm[,-c(3,9)]
head(red_final_mrna)
rownames(red_final_mrna) <- red_final_mrna$GeneSymbol
red_final_mrna$GeneSymbol <- NULL
head(red_final_mrna)
red_final_mrna <- red_final_mrna[,c(2,3,8,5,6,7,1,4)]

design <- data.frame(row.names= colnames(red_final_mirna),sample= colnames(red_final_mirna),condition=c(rep("Responder",3),rep("NonResponder",3)))
design$condition <- factor(design$condition,levels=c("NonResponder","Responder"))
design.matrix <- model.matrix(~ condition,data=design)
keepNonxR2 <- which((rowSums(red_final_mrna[,c(1:3)] > 10) >= 2 ) | (rowSums(red_final_mrna[,c(4:6)] > 10) >= 2))
head(red_final_mirna_filt)
red_final_mrna_filt <- red_final_mrna[names(keepNonxR2),]
length(keepNonxR)
y_mrna <- DGEList(counts=red_final_mrna_filt,genes=row.names(red_final_mrna_filt))
y_mrna <- estimateCommonDisp(y_mrna,method="TMM")
plotMDS(y_mrna)
y_mrna2 <- y_mrna[,-c(7,8)]
y_mrna2<- estimateGLMTrendedDisp(y_mrna2,design.matrix)
y_mrna2 <- estimateGLMTagwiseDisp(y_mrna2,design.matrix)
fit_mrna <- glmFit(y_mrna2,design.matrix)
resmrna <- glmLRT(fit_mrna,coef=2)
head(resmirna)
results_mrna <- topTags((resmrna),n=nrow(y_mrna2$rows))
summary(results_mrna$table$FDR)
rownamesPvaluem <- rownames(subset(results_mrna$table,results_mirna$table$PValue < 0.05))
head(rownamesPvalue)
y_mrna3 <- y_mrna[rownamesPvaluem,]
plotMDS(y_mrna3)

cpmymrna <- cpm(y_mrna3)
heatmap.2(log2(cpmymrna+1/(rowMeans(cpmymrna) +1)))
dim(cpmyn)
library(gplots)
pdf(heatmap.2(cpmyn),file=heatmap.pdf)
countsmrna <- cbind(y_mrna2$counts,cpm(y_mrna2))
resultsmrna <- cbind(countsmrna,results_mrna$table[rownames(countsmrna),])
write.table (resultsmrnac, file="excerpt_mrna_NOTcorrected.tsv",col.names = T,row.names=T,sep="\t")
#####correction
groupsl <- c(rep(1,3),rep(2,5))
batchesl <- c(1,2,2,2,1,1,1,2)
correctmr <- ComBat_seq(counts = as.matrix(y_mrna$counts),batch=batchesl,group=groupsl)
head(correctl)
y_mrnac <- DGEList(counts=correctmr,genes=row.names(correctmr))
dds <- DESeqDataSetFromMatrix(countData = round(correctmr),
                              colData = design2,
                              design =~ condition)

dds <- DESeq(dds)
vsdmr <- vst(dds)

plotMDS(y_mrnac)
y_mrnac <- estimateCommonDisp((y_mrnac))
cpm_mrnac <- log2(cpm(y_mrnac)+1)
dist_mrnac <- dist(t(cpm_mrnac), method = 'euclidean')
hclust_mrnac <- hclust(pmds$distance.matrix.squared)
head(pmds$distance.matrix.squared)
head(dist_mrnac)
plot(hclust_mrnac)
cmdscl <- cmdscale(dist_mrnac, eig = TRUE, k = 2)
pmds <- plotMDS(y_mrnac)
head(pmds)
xcl <- pmds$eigen.vectors[,1]
ycl <- pmds$eigen.vectors[, 2]
#yl2
colnames(yl2$count)
head(cmdscl$eig)
head(pmds)
#yl
ggplot() + geom_point(data = as.data.frame(pmds$eigen.values) , mapping = aes(x = xcl, y = ycl), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mrnac) ,cmdscl$eig), mapping = aes(x=xcl, y=ycl , label = colnames(cpm_mrnac)  )) + labs(title = "MDS_lung")
#design2mi <- data.frame(row.names= colnames(correctm),sample= colnames(correctm),condition=c(rep("Responder",3),rep("NonResponder",5)))
#design2$condition <- factor(design2$condition,levels=c("NonResponder","Responder"))
#design.matrix2 <- model.matrix(~ condition,data=design2)
y_mrnac<- estimateGLMTrendedDisp(y_mrnac,design.matrix2)
y_mrnac <- estimateGLMTagwiseDisp(y_mrnac,design.matrix2)
fit_mrnac <- glmFit(y_mrnac,design.matrix2)
resmrnac <- glmLRT(fit_mrnac,coef=2)
head(resmirna)
library(ggplot2)
library(ggrepel)
results_mrnac <- topTags((resmrnac),n=nrow(y_mrnac$rows))
dim(subset(results_mrnac$table, results_mrnac$table$FDR < 0.05))
rownamesFDRmrnac <- rownames(subset(results_mrnac$table,results_mrnac$table$FDR < 0.05))
head(rownamesPvalue)
y_mrnac3 <- y_mrnac[rownamesFDRmrnac,]
plotMDS(y_mrnac3)
pmds2 <- plotMDS(y_mrnac3)
head(pmds)
xcl2 <- pmds2$eigen.vectors[,1]
ycl2 <- pmds2$eigen.vectors[, 2]
ggplot() + geom_point(data = as.data.frame(pmds2$eigen.values) , mapping = aes(x = xcl2, y = ycl2), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mrnac) ,pmds2$eigen.values), mapping = aes(x=xcl2, y=ycl2 , label = colnames(cpm_mrnac)  )) + labs(title = "MDS_lung")
countsmrnac <- cbind(y_mrnac$counts,cpm(y_mrnac))
head(countsl)
head(totall)
resultsmrnac <- cbind(countsmrna,results_mrnac$table[rownames(countsmrnac),])
head(resultsmrnac)
library(EnhancedVolcano)
colnames(resultsmrnac)
resultsmrnacV <- resultsmrnac[,17:22]
resultsmrnacV$diffexpressed <- "NO"
resultsmrnacV$diffexpressed[resultsmrnacV$logFC > 2 & resultsmrnacV$FDR < 0.05] <- "UP"
resultsmrnacV$diffexpressed[resultsmrnacV$logFC < -2 & resultsmrnacV$FDR < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
ggplot(data=resultsmrnacV, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) 
fdrmrna <- subset(resultsmrnac[,9:16],resultsmrnac$FDR < 0.05)
colnames(fdrmrna)
heatmap.2(as.matrix(log2((fdrmrna+1)/rowMeans(fdrmrna+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="mRNA Excerpt DE")
write.table (resultsmrnac, file="excerpt_mrna_corrected.tsv",col.names = T,row.names=T,sep="\t")

##################################excerpt mirna small rna seq
mirna_list2 <- list.files(path="../smallExceprt/",recursive = T, pattern="miRNAmature_sense.txt")
mirna_list2
finalmir <- read.table("../smallExceprt/100-1-RNA-Lung-50ng-NovaSeq_S11_L001_R1_001.fastq/readCounts_miRNAmature_sense.txt",header=T)
head(finalmir)
finalmir <- finalmir[,-c(2,3,5)]
head(temp)
fpath
for (i in 2:length(mirna_list2)) {
  fpath <- paste("../smallExceprt/",mirna_list2[i],sep="")
  temp <- read.table(fpath,header=T)
  temp <- temp[,-c(2,3,5)]
  finalmir <- merge(finalmir,temp,by="ReferenceID",all=T)
}
head(finalmir)  
namesn <- colnames(tableall)[4:13]
namesn <- namesn[-7]
namesn[1] <- gsub("_100ng","",namesn[1])
mirna_list2

colnames(finalmir)[2:9] <- namesn

finalmir[is.na(finalmir)] <- 0

head(red_final_mirna)
rownames(finalmir) <- finalmir$ReferenceID
finalmir$ReferenceID <- NULL
finamir <- finalmir[,c(2,3,8,5,6,7,1,4)]
design2
head(finamir)

keepNonxR <- which((rowSums(finamir[,c(1:3)] > 10) >= 2 ) | (rowSums(finamir[,c(4:6)] > 10) >= 2))
head(red_final_mirna_filt)
finamir <- finamir[names(keepNonxR),]
head(finamir)

correctmi <- ComBat_seq(counts = as.matrix(finamir),batch=batchesl,group=groupsl)
head(correctmi)
y_mirnaSc <- DGEList(counts=correctmi,genes=row.names(correctmi))
plotMDS(y_mirnaSc)
y_mirnaS <- DGEList(counts=finamir,genes=row.names(finamir))
plotMDS(y_mirnaS)
#design2mi <- data.frame(row.names= colnames(correctm),sample= colnames(correctm),condition=c(rep("Responder",3),rep("NonResponder",5)))
#design2$condition <- factor(design2$condition,levels=c("NonResponder","Responder"))
#design.matrix2 <- model.matrix(~ condition,data=design2)
y_mirnaSc<- estimateGLMTrendedDisp(y_mirnaSc,design.matrix2)
y_mirnaSc <- estimateGLMTagwiseDisp(y_mirnaSc,design.matrix2)
fit_mirnaSc <- glmFit(y_mirnaSc,design.matrix2)
resmirnaSc <- glmLRT(fit_mirnaSc,coef=2)
head(resmirna)

results_mirnaSc <- topTags((resmirnaSc),n=nrow(y_mirnaSc$rows))
rownamesFDRmirnaSc <- rownames(subset(results_mirnaSc$table,results_mirnaSc$table$PValue < 0.05))
head(rownamesFDRmirnaSc)
y_mirnaSc3 <- y_mirnaSc[rownamesFDRmirnaSc,]
y_mirnaSc3
plotMDS(y_mirnac3)
pmdSs3 <- plotMDS(y_mirnaSc3)
head(pmds)
xcl5 <- pmdSs3$eigen.vectors[,1]
ycl5 <- pmdSs3$eigen.vectors[, 2]
cpm_mirnaSc <- cpm(y_mirnaSc)
ggplot() + geom_point(data = as.data.frame(pmdSs3$eigen.values) , mapping = aes(x = xcl5, y = ycl5), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mirnaSc) ,pmdSs3$eigen.values), mapping = aes(x=xcl5, y=ycl5 , label = colnames(cpm_mirnaSc)  )) + labs(title = "MDS_mirna_lung")
countsmirnaS <- cbind(y_mirnaSc$counts,cpm(y_mirnaSc))
head(countsl)
head(totall)
resultsmirnaSc <- cbind(countsmirnaS,results_mirnaSc$table[rownames(countsmirnaS),])
head(resultsmrnac)
library(EnhancedVolcano)

resultsmirnaScV <- resultsmirnaSc[,17:22]
resultsmirnaScV$diffexpressed <- "NO"
resultsmirnaScV$diffexpressed[resultsmirnaScV$logFC > 2 & resultsmirnaScV$PValue < 0.05] <- "UP"
resultsmirnaScV$diffexpressed[resultsmirnaScV$logFC < -2 & resultsmirnaScV$PValue < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
ggplot(data=resultsmirnaScV, aes(x=logFC, y=-log10(PValue), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) 
fdrmirnaS <- subset(resultsmirnaSc[,9:16],resultsmirnaSc$PValue < 0.05)
colnames(fdrmrna)
rownames(fdrmirna)

rownames(fdrmirna) <- gsub(":MIMAT\\d+:Homo:sapiens:miR.+","",rownames(fdrmirna))
heatmap.2(as.matrix(log2((fdrmirna+1)/rowMeans(fdrmirna+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,12), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="miRNA Excerpt corrected DE")
write.table (resultsmirnaSc, file="excerpt_mirna_SMALLRNA_corrected.tsv",col.names = T,row.names=T,sep="\t")



barplot(y_der$samples$lib.size,names.arg=rownames(y_der$samples))
plot(log2(count100red))
plotMDS(y_der2)
pdf("report.pdf")
barplot(y_der$samples$lib.size,names.arg=rownames(y_der$samples),main="derfinder")
plotMDS(y_der2, main="derfinder only known samples before correction")
plotMDS(y_der3,main="derfinder all samples with fdr sig NonxResponders")
ggplot() + geom_point(data = as.data.frame(pmds2$eigen.values) , mapping = aes(x = xcl2, y = ycl2), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mrnac) ,pmds2$eigen.values), mapping = aes(x=xcl2, y=ycl2 , label = colnames(cpm_mrnac)  )) + labs(title = "MDS Derfinder Corrected")
ggplot(data=resultsdercV, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  labs(title = "Volcano Derfinder Corrected")+
  scale_color_manual(values=c("blue", "black", "red"))
  
heatmap.2(as.matrix(log2((fdrderc+1)/rowMeans(fdrderc+1))),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),labRow=F,main="Derfinder corrected DE")
ggplot() + geom_point(data = as.data.frame(pmds2$eigen.values) , mapping = aes(x = xcl2, y = ycl2), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mrnac) ,pmds2$eigen.values), mapping = aes(x=xcl2, y=ycl2 , label = colnames(cpm_mrnac)  )) + labs(title = "MDS miRNA Excerpt Corrected")
ggplot(data=resultsmirnacV, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  labs(title = "Volcano miRNA Excerpt Corrected")+
  scale_color_manual(values=c("blue", "black", "red")) 
heatmap.2(as.matrix(log2((fdrmrna+1)/rowMeans(fdrmrna+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="mRNA Excerpt DE")
ggplot() + geom_point(data = as.data.frame(pmds2$eigen.values) , mapping = aes(x = xcl2, y = ycl2), color = c(rep("green",3),rep("red",3),"red","red"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpm_mrnac) ,pmds2$eigen.values), mapping = aes(x=xcl2, y=ycl2 , label = colnames(cpm_mrnac)  )) + labs(title = "MDS mRNA Excerpt Corrected")
ggplot(data=resultsmrnacV, aes(x=logFC, y=-log10(FDR), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  labs(title = "Volcano mRNA Excerpt Corrected")+
  scale_color_manual(values=c("blue", "black", "red")) 

heatmap.2(as.matrix(log2((fdrmrna+1)/rowMeans(fdrmrna+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",5)),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="mRNA Excerpt DE")
dev.off()
