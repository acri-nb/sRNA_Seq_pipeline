tableall <- read.table("lung.Output.txt")
colnames(tableall) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(tableall))
rownames(tableall) <- paste(tableall$contig,":",tableall$start,"-",tableall$end,sep="")
countl <- round(as.data.frame(tableall[,4:13]))
dim(countl)
head(spikeinl)
rownames(spikeinl) <- spikeinl$spikein
spikeinl$spikein <- NULL
head(spikein[,colnames(count)])
head(count)
count2l <- rbind(countl,spikeinl[,colnames(countl)])
colnames(countl) <- c("BTD100","BTD122","BTD131","BTD193","BTD198","BTD36","ACRI500","BTD75","BTD83","BTD97")

head(countl)

count_red <- countl[names(which((rowSums(countl[,c(5,6,7,9)] > 10) >=2) | rowSums(countl[,c(2,3,4,8,10)] > 10) >=3 )),]  
head(count_red)
yl <- DGEList(counts=count_red[,-c(1,4,6)],genes=row.names(count_red[,-c(1,4,6)]))
colnames(count_red)
head(y2l$counts)
cor(spikeinl)
cor(cpm(yl))
keep <- rowSums(cpm(yl)>1) >= 8
y2l <- yl[keep,]
yl <- calcNormFactors(yl, method="TMM")
yl <- estimateCommonDisp(yl)
plotMDS(yl)
rownames(y2l)
spikein <- cpm(y2l)[959:1000,]
rest <- cpm(y2l)[1:958,]
plotMDS(yl)
colnames(yl$counts)
cpms <- cpm(yl)
head(cpms)
dist_matl2 <- dist(t(cpms[,-7]), method = 'euclidean')
hclust_avgl2 <- hclust(dist_matl2, method = 'complete')
pdf("lungs_two.pdf")
plot(hclust_avgl2)
cmdsl2 <- cmdscale(dist_matl2, eig = TRUE, k = 2)
xl2 <- cmdsl2$points[, 1]
yl2 <- cmdsl2$points[, 2]
ggplot() + geom_point(data = as.data.frame(cmdsl2$eig) , mapping = aes(x = xl2, y = yl2), color = c(rep("green",3),rep("red",3),"red","green"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpms)[-7] ,cmdsl2$eig), mapping = aes(x=xl2, y=yl2 , label = colnames(cpms)[-7]  )) + labs(title = "Lung MDS")
dev.off()
pdf("lung_clustering.pdf")
dev.off()
groupsl <- c(rep(1,3),rep(2,3),1,2,1)
batchesl <- c(1,2,2,2,1,2,2,1,1)
groupsl <- c(rep(1,3),rep(2,3),1,2,1)
batchesl <- c(1,2,2,2,1,2,2,1,1)
correctl <- ComBat_seq(counts = as.matrix(count_red[,-c(1)]),batch=batchesl,group=groupsl)
plotMDS(ycl)
ycl <- DGEList(counts=correctl,genes=row.names(correctl))
barplot(ycl$samples$lib.size)
ycl <- calcNormFactors(ycl, method="TMM")
cpmsc <- cpm(ycl)
dist_matcl2 <- dist(t(cpmsc), method = 'euclidean')
hclust_avgcl2 <- hclust(dist_matcl2, method = 'complete')
pdf("corrected.pdf")
plot(hclust_avgcl2)
dev.off()
cmdscl <- cmdscale(dist_matcl2, eig = TRUE, k = 2)
xcl <- cmdscl$points[, 1]
ycl <- cmdscl$points[, 2]
#yl2
colnames(yl2$count)
#yl
ggplot() + geom_point(data = as.data.frame(cmdscl$eig) , mapping = aes(x = xcl, y = ycl), color = c(rep("green",3),rep("red",3),"green","red","green"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(cpmsc) ,cmdscl$eig), mapping = aes(x=xcl, y=ycl , label = colnames(cpmsc)  )) + labs(title = "MDS_lung")
dev.off()



designl <- data.frame(row.names= colnames(yl$counts),sample= colnames(yl$counts),condition=c(rep("respondent",3),rep("nonRespondent",3),"respondent","nonRespondent","respondent"))
designl$condition <- factor(designl$condition,levels=c("respondent","nonRespondent"))
designl
design.matrixl <- model.matrix(~ condition,data=designl)
yl
design.matrixl
yl2 <- estimateGLMCommonDisp(yl,design.matrixl)
nrow(designl)
ncol(yl)
design
barplot(yl$samples$lib.size)
yl2 <- estimateGLMTrendedDisp(yl2,design.matrixl)
yl2 <- estimateGLMTagwiseDisp(yl2,design.matrixl)
fitl <- glmFit(yl2,design.matrixl)
resl <- glmLRT(fitl,coef=2)
topTags((resl))
yl2["chrX:19376772-19376804",]

totall <- topTags(resl,n = nrow(yl2$counts))
dim(subset(totall$table, totall$table$PValue < 0.05))
countsl <- cbind(yl2$counts,cpm(yl2))
head(countsl)
head(totall)
resultsl <- cbind(countsl,totall$table[rownames(countsl),])
fdr <- subset(resultsl[,10:18],resultsl$FDR < 0.05)
fdr2 <- subset(resultsl[,10:18],resultsl$PValue < 0.05)
pdf ( "volcano.pdf")

with(resultsl,plot(logFC, -log10(PValue),pch=20, main="Volcano plot"))
with(subset(resultsl,PValue<.05),points(logFC, -log10(PValue),pch=20,col="purple"))
with(subset(resultsl,logFC > 1.5 | logFC < -1.5), points(logFC, -log10(PValue),pch=20, col="orange"))
with (subset(resultsl,PValue < 0.05 & (logFC > 1.5 | logFC < -1.5)), points(logFC, -log10(PValue), pch=20,col="yellow"))
with (resultsl[list,], text(logFC, -log10(PValue),genes))
resultsl$status = NA
head(resultsl$status)
resultsl[list3,]$status <- tableall[list3,]$ids
head(list3)
#resultsl[list,]$status=TRUE
library(ggrepel)
resultsl$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
resultsl$diffexpressed[resultsl$logFC > 1.5 & resultsl$PValue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
resultsl$logFC < -1.5
resultsl$logFC < -1.5 & resultsl$PValue < 0.05
resultsl$diffexpressed[resultsl$logFC < -1.5 & resultsl$PValue < 0.05] <- "DOWN"
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
table(resultsl$diffexpressed)
colnames(resultsl)
resultsld <- resultsl [!duplicated(resultsl[c(20,23)]),]
resultsld[!duplicated(resultsld$logFC),]
resultsd2 <- resultsl[,20:26]
table(resultsld$status)
ggplot(data=resultsd2, aes(x=logFC, y=-log10(PValue), col=diffexpressed, label=status)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size=4,box.padding = unit(0.5, "lines")) +
  scale_color_manual(values=c("blue", "black", "red")) 
#  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
 # geom_hline(yintercept=-log10(0.05), col="red")
pdf("volcanoplotsv3.pdf",width=10)
EnhancedVolcano(resultsd2,
                lab = resultsd2$status,
                x = 'logFC',
                y = 'PValue',
                ylim = c(0,5.8),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 5e-2,
                FCcutoff = 1.5,
                pointSize = 2.0,
                labSize = 2.5,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
         #       boxedLabels = TRUE,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.80)
dev.off()
head(resultsl)
list=c("chr5:181241818-181241889","chrX:19376772-19376804","chr15:25051477-25051569","chr2:148881724-148881821","chr20:2654211-2654286","chr19:53671977-53671999")
list2 =c("chr7:10982594-10982669","chr12:79795723-79795764","chr1:155679145-155679240","chr10:17437504-17437522","chr3:186786673-186786767")
list3 = c(list,list2)
tableall[list,]$ids
dev.off()
dist_matl <- dist(t(fdr), method = 'euclidean')
cmdsl <- cmdscale(dist_matl, eig = TRUE, k = 2)
xl <- cmdsl$points[, 1]
yl <- cmdsl$points[, 2]
#yl2
colnames(yl2$count)
#yl
ggplot() + geom_point(data = as.data.frame(cmdsl$eig) , mapping = aes(x = xl, y = yl), color = c(rep("green",3),rep("red",3),"green","red","green"), alpha = 0.5) + geom_text_repel(data =data.frame(colnames(yl2$counts) ,cmdsl$eig), mapping = aes(x=xl, y=yl , label = colnames(yl2$counts)  )) + labs(title = "figure5:  MDS configuration of Japan Prefectures with labels")
cmdsl                                                                                                                                                                                        
colnames(tableall)
head(resultsl)
resultsl2 <- cbind(resultsl,tableall[rownames(resultsl),14:19])
nocorrected_pvalue_sub <- resultsl2[which(resultsl2$PValue < 0.05),]
protein_notc<- nocorrected_pvalue_sub [which(grepl("protein_coding",nocorrected_pvalue_sub $type)),]
write.table (protein_notc, file="protein_notcorrected.tsv",col.names = T,row.names=T,sep="\t")
write.table(resultsl2,"lung_comparision.tsv",sep="\t")
library(gplots)
pheatmap(log2(fdr + 1))
heatmap.2(as.matrix(log2((fdr+1)/rowMeans(fdr+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",3),"green","red","green"),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="No corretion by FDR < 0.05")
heatmap.2(as.matrix(log2((fdr2+1)/rowMeans(fdr2+1))+0.5),trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(59),ColSideColors = c(rep("green",3),rep("red",3),"green","red","green"),margins=c(8,8), breaks=c(seq(-2,-1,length=20),seq(-0.9,0.9,length=20), seq(1,3,length=20)),main="No corretion by PVALUE < 0.05")
head(fdr2)
ycl
ycl2 <- estimateGLMCommonDisp(ycl,design.matrixl)
nrow(designl)
ncol(yl)
design
barplot(yl$samples$lib.size)
ycl2 <- estimateGLMTrendedDisp(ycl2,design.matrixl)
ycl2 <- estimateGLMTagwiseDisp(ycl2,design.matrixl)
fitcl <- glmFit(ycl2,design.matrixl)
rescl <- glmLRT(fitcl,coef=2)
resultsc <- topTags((rescl),n=nrow(ycl2$rows))
topTags(rescl)
fdrc <- subset(resultsc$table,resultsc$table$FDR <  0.05)
fdrc2 <- subset(resultsc$table,resultsc$table$PValue <  0.05)
head(fdrc)
cpmfdr <- as.matrix(cpm(ycl2)[rownames(fdrc),])
cpmfdr2 <- as.matrix(cpm(ycl2)[rownames(fdrc2),])
library(gplots)
cpmfdr
pdf("heatmaps_lung.pdf")

heatmap.2(log2((cpmfdr+1)/rowMeans(cpmfdr+1))+0.5,trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(249),ColSideColors = c(rep("green",3),rep("red",3),"green","red","green"),margins=c(8,8), breaks=c(seq(-2,-1,length=50),seq(-0.9,0.9,length=150), seq(1,2,length=50)),main="with correction and fdr < 0.05")
heatmap.2(log2((cpmfdr+1)/rowMeans(cpmfdr+1))+0.5,trace='none',symm=F,symkey=F,symbreaks=T,col = bluered(249),ColSideColors = c(rep("green",3),rep("red",3),"green","red","green"),margins=c(8,8), breaks=c(seq(-2,-1,length=50),seq(-0.9,0.9,length=150), seq(1,2,length=50)),main="with correction and pvalue < 0.05")
dev.off()
resultstc2 <- cbind(ycl2$counts,cpm(ycl2))
resultstc2 <- cbind(resultstc2,resultsc$table[rownames(resultstc2),])
colnames(tableall)
resultstc2 <- cbind(resultstc2,tableall[rownames(resultstc2),14:19])
colnames(resultstc2)
corrected_pvalue_sub <- resultstc2[which(resultstc2$PValue < 0.05),]
protein_corrected <- corrected_pvalue_sub[which(grepl("protein_coding",corrected_pvalue_sub$type)),]
write.table (protein_corrected, file="protein_corrected.tsv",col.names = T,row.names=T,sep="\t")
write.table (resultstc2, file="results_GBM_groups.tsv",col.names = T,row.names=T,sep="\t")
designl
head(resultstG[which(resultstG$FDR < 0.05),])
reactome <- read.csv2("~/Downloads/Panther_2016_table (2).txt",sep="\t")
head(reactome)
splits <- strsplit(reactome$Overlap, "/")
splits[[1]][1]
length(splits)
numbers <- ""
for (i in 1:length(splits)){
 numbers <-c (numbers,(splits[[i]][1])) 
  
}
head(numbers)
numbers <- numbers[-1]
reactome$size <- as.numeric(numbers)
head(reactome)
sub_reactome <- subset(reactome,reactome$P.value < 0.1)
plot(sub_reactome$size,type="bar")
barplot(sub_reactome$size,horiz=T)
head(sub_reactome$size)
table_reactome <- data.frame(row.names=sub_reactome$Term,size=sub_reactome$size
)    
str(sub_reactome)
barplot(t(table_reactome[c(26:1)],),horiz=T)
sub_reactome$Term <- factor(sub_reactome$Term,levels=rev(sub_reactome$Term))
ggplot(sub_reactome, aes(x = Term, y = size)) +   geom_col(size = 1, color = "darkblue", fill = "white") + coord_flip()

sub_reactome$P.value <- as.numeric(sub_reactome$P.value) 
ggplot(sub_reactome,aes(x = Term, y = P.value)) + geom_line(size = 1.5, color="red", group = 1) + coord_flip()
-log10(sub_reactome$P.value)
pdf("third_version_barplot.pdf")
barplot(sub_reactome$size)
ggplot(sub_reactome) + 
  geom_col(aes(x = Term, y = size), size = 0.5, color = "black", fill = "red") +
  geom_line(aes(x = Term, y = -log10(P.value)), size = 0.5, color="black", group = 1) +  
  geom_point(aes(x = Term, y = -log10(P.value)), size = 1, color="black", group = 1) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "-log10 P.value")) + 
  ylab("Number of genes") +
  coord_flip()+ theme_bw()
dev.off()
sub_reactome$P.value
-log10(0.05)
