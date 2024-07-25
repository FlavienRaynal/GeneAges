## Chromatin network rearrangements in differentiation and cancer reveal connections between evolutionary processes and gene regulation


##############
##  Fig 1A  ## # Boxplot Trigos ages / Monocyte Mean Exp
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(ggplot2)
library(ggpubr)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
ev=read.table(paste0(datasets_dir,"Exp_EV_Ecker.txt"),header = T) # Exp var

### Merge datasets
trigos=merge(trigos,ev,by="gene_id")

### Make Boxplot
color_mono="#B3CDE3"
rm(p)
p=ggplot(trigos, aes(x = Phylostrata, y = mean_mono, group=Phylostrata)) 
p= p + geom_boxplot(coef=6,fill = color_mono,lwd=1.1)
p= p + ylab("Monocyte Mean Exp")
p= p + theme_classic() 
p= p + stat_summary(fun = mean, colour = "red", geom="line")
p= p + stat_summary(aes(x=Phylostrata,y=mean_mono), fun=mean, colour="red", geom="line",group=1,linetype = 2,size=1.2)
p= p + stat_summary(aes(x=Phylostrata,y=mean_mono), fun=mean, colour="red", geom="point",group=1)
print(p)






##############
##  Fig 1B  ## # Boxplot Trigos ages / EV in Monocyte
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(ggplot2)
library(ggpubr)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
ev=read.table(paste0(datasets_dir,"Exp_EV_Ecker.txt"),header = T) # Exp var

### Merge datasets
trigos=merge(trigos,ev,by="gene_id")
trigos=trigos[which(trigos$ev_mono>(-10)),] #remove outliers

### Make Boxplot
color_mono="#B3CDE3"
rm(p)
p=ggplot(trigos, aes(x = Phylostrata, y = ev_mono, group=Phylostrata)) 
p= p + geom_boxplot(coef=6,fill = color_mono,lwd=1.1)
p= p + ylab("Monocyte EV") #for mean
p= p + theme_classic() 
p= p + stat_summary(fun = mean, colour = "red", geom="line")
p= p + stat_summary(aes(x=Phylostrata,y=ev_mono), fun=mean, colour="red", geom="line",group=1,linetype = 2,size=1.2)
p <- p + stat_summary(aes(x=Phylostrata,y=ev_mono), fun=mean, colour="red", geom="point",group=1)
print(p)






##############
##  Fig 1E  ## # Polycomb target genes in Monocyte
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
library(genomation)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
k27me3=readBed(paste0(datasets_dir,"H3K27me3_Monocyte_hg19_peaks.narrowPeak")) # Gene Ages

### Get age-annotated gene promoters coordinates (TSS +/- 1,5 kb)
genes=genes(EnsDb.Hsapiens.v75)
seqlevels(genes)=paste0("chr",seqlevels(genes))
genes=genes[which(genes$gene_name %in% trigos$GeneID)]
merge=unique(merge(as.data.frame(trigos)[,c(1,4,5)],as.data.frame(genes)[,c(1:3,5,7)],by.x="GeneID",by.y="gene_name"))
merge=GRanges(merge$seqnames,IRanges(merge$start,merge$end),strand=merge$strand,GeneID=merge$GeneID,Phylostrata=merge$Phylostrata,mainAge=merge$mainAge)
prom=merge
start=ifelse(strand(prom)=="+",start(prom)-1500,end(prom)-1500)
end=ifelse(strand(prom)=="+",start(prom)+1500,end(prom)+1500)
start(prom)=start;end(prom)=end

### Get targeted promoters
ol=findOverlaps(prom,k27me3)
target_prom=as.data.frame(prom[queryHits(ol)],row.names = NULL)
target_prom=unique(target_prom)

### Compute percentages
tab_3ages=table(target_prom$mainAge)[c(3,1,2)]
tab_3ages=rbind(tab_3ages,table(trigos$mainAge)[c(3,1,2)])
tab_3ages[1,]=tab_3ages[1,]/tab_3ages[2,]*100
tab_3ages[2,]=50-(tab_3ages[1,])

### Make Barplot
col=c("red3","grey")
b=barplot(tab_3ages,col=col, ylab="",main="Monocyte Polycomb target genes \n Trigos Gene Ages - Prom",ylim=c(0,50),xaxt="n",axes=0,border = "white") #3mainAges
text(b, y= tab_3ages[1,], labels = paste0(round(tab_3ages[1,]),"%"),cex=1.6,pos=3)
mtext("Polycomb Targets (%)", side = 2, line = 2.5,cex=1.4)
axis(2, at=c(0,25,50),las=2,cex.axis=1.2,lwd=1.2,font=1)
text(cex=1.6, x=b+.4, y=-7, colnames(tab_3ages), xpd=T,pos=2)





##############
##  Fig 1F  ## # Average number of CpGs by gene age
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
library(genomation)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
#meth=read.table(paste0(datasets_dir,"Meth_scores_Monocyte_Chen.txt"),header = T) # Methylation scores
meth=readRDS("/Users/fla/Documents/Work/GeneAges/EV/monocytes-meth-feature-table.rds")

### Get age-annotated gene promoters coordinates (from TSS-2kb to TSS+200bp)
genes=genes(EnsDb.Hsapiens.v75)
prom=genes[which(genes$gene_name %in% trigos$GeneID),]
seqlevels(prom)=paste0("chr",seqlevels(prom))
start=ifelse(strand(prom)=="+",start(prom)-2000,end(prom)-200)
end=ifelse(strand(prom)=="+",start(prom)+200,end(prom)+2000)
start(prom)=start;end(prom)=end

### Get genes overlapping with CpG
cpg=GRanges(meth$chrom,IRanges(meth$start,meth$end))
ol=findOverlaps(cpg,prom)
meth$gene=NA
meth$gene[queryHits(ol)]=prom$gene_name[subjectHits(ol)] #annotate cpg overlapping genes

### Get number of CpG by gene
nb=data.frame(tab=table(meth$gene),"NA"=NA)
nb=merge(nb,trigos[,c(1,4,5)],by.x="tab.Var1",by.y="GeneID") #add gene ages
nb=unique(nb)

### Average for 3 main ages
m=vector()
for(i in c("UC","EM","MM")){
  m=c(m,mean(nb[which(nb$mainAge==i),2]))
}
names(m)=c("UC","EM","MM")

### Make Barplot
color_3ages=c("#f7766d","#02bb38","#609bff")  
b=barplot(m,col=color_3ages, ylab="",main="Monocyte Number of CpG on Promoters\n Trigos Gene Ages",ylim=c(0,8.3),xaxt="n",axes=0,border = "white")
pos=ifelse(m>0,round(m,1),"")
text(b, m, labels=pos, cex=1.6,pos=3)
mtext("", side = 2, line = 2,cex=1.4)
axis(2, at=c(0,8),las=2,cex.axis=1.6,lwd=1.2,font=1)
text(x=b, y=-1, names(m), cex=1.6,xpd=T)
cor.test(1:3,m,method="spearman") #Sperman corr





##############
##  Fig 1G  ## # Average methylation scores
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
#meth=read.table(paste0(datasets_dir,"Meth_scores_Monocyte_Chen.txt"),header = T) # Methylation scores
meth=readRDS("/Users/fla/Documents/Work/GeneAges/EV/monocytes-meth-feature-table.rds")

### Get age-annotated gene promoters coordinates (from TSS-2kb to TSS+200bp)
genes=genes(EnsDb.Hsapiens.v75)
prom=genes[which(genes$gene_name %in% trigos$GeneID),]
seqlevels(prom)=paste0("chr",seqlevels(prom))
start=ifelse(strand(prom)=="+",start(prom)-2000,end(prom)-200)
end=ifelse(strand(prom)=="+",start(prom)+200,end(prom)+2000)
start(prom)=start;end(prom)=end

### Get genes overlapping with CpG
meth$mean=rowMeans(meth[,4:ncol(meth)])
cpg=GRanges(meth$chrom,IRanges(meth$start,meth$end))
ol=findOverlaps(cpg,prom)
meth$gene=NA
meth$gene[queryHits(ol)]=prom$gene_name[subjectHits(ol)] #annotate cpg overlapping genes

### Get number of CpG by gene
nb=data.frame(tab=table(meth$gene),"NA"=NA)
nb=merge(nb,trigos[,c(1,4,5)],by.x="tab.Var1",by.y="GeneID") #add gene ages
nb=unique(nb)

### Get DNA methylation frequency 
freq=aggregate(. ~ gene, data = meth[,c(200,201)], mean)
freq=merge(freq,trigos[,c(1,4,5)],by.x="gene",by.y="GeneID")
freq=unique(freq)

### Average for 16 ages (supp?)
# m=vector()
# for(i in 1:16){
#   m=c(m,mean(freq[which(freq$Phylostrata==i),2]))
# }
# names(m)=paste0("Class",1:16)

### Make Boxplot
rm(p)
freq$mainAge=factor(nb$mainAge,levels = c("UC","EM","MM"))
color_3ages=c("#f7766d","#02bb38","#609bff")  
p=ggplot(freq, aes(x = mainAge, y = mean, group=mainAge)) 
p= p + geom_boxplot(coef=6,fill = color_3ages,lwd=1.1)
p= p + ylab("Methylation Scores") #for mean
p= p + ggtitle("Methylation Scores - 3 main ages") #for mean
p= p + theme_classic() 
p= p + stat_summary(fun = mean, colour = "red", geom="line")
p= p + stat_summary(aes(x=mainAge,y=mean), fun=mean, colour="red", geom="line",group=1,linetype = 2,size=1.2)
p <- p + stat_summary(aes(x=mainAge,y=mean), fun=mean, colour="red", geom="point",group=1)
print(p)





##############
##  Fig 1H  ## # Average DNA methylation variability
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
library(ggplot2)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
meth=readRDS("/Users/fla/Documents/Work/GeneAges/EV/monocytes-meth-feature-table.rds")
mv=read.table(paste0(datasets_dir,"Meth_Var_Monocyte.txt"),header = T) # Methylation Variability

### Get age-annotated gene promoters coordinates (from TSS-2kb to TSS+200bp)
genes=genes(EnsDb.Hsapiens.v75)
prom=genes[which(genes$gene_name %in% trigos$GeneID),]
seqlevels(prom)=paste0("chr",seqlevels(prom))
start=ifelse(strand(prom)=="+",start(prom)-2000,end(prom)-200)
end=ifelse(strand(prom)=="+",start(prom)+200,end(prom)+2000)
start(prom)=start;end(prom)=end

### Get genes overlapping with CpG
meth$mean=rowMeans(meth[,4:ncol(meth)])
cpg=GRanges(meth$chrom,IRanges(meth$start,meth$end))
ol=findOverlaps(cpg,prom)
meth$gene=NA
meth$gene[queryHits(ol)]=prom$gene_name[subjectHits(ol)] #annotate cpg overlapping genes

### Annotate with meth var and gene ages
meth$mv=mv$mv
meth_mv=aggregate(. ~ gene, data = meth[,c(201,202)], mean)
meth_mv=merge(meth_mv,trigos[,c(1,4,5)],by.x="gene",by.y="GeneID")
meth_mv=unique(meth_mv)

### Make Boxplot
meth_mv$mainAge=factor(meth_mv$mainAge,levels = c("UC","EM","MM"))
color_3ages=c("#f7766d","#02bb38","#609bff")  
rm(p)
p=ggplot(meth_mv, aes(x = mainAge, y = mv, group=mainAge)) 
p= p + geom_boxplot(coef=6,fill = color_3ages,lwd=1.1,outlier.shape = NA)
p= p + ylab("Methylation Variability") 
p= p + ggtitle("Methylation Variability - 3 main ages") #for mean
p= p + theme_classic() + ylim(-1.5,1.5)
p= p + stat_summary(fun = mean, colour = "red", geom="line")
p= p + stat_summary(aes(x=mainAge,y=mv), fun=mean, colour="red", geom="line",group=1,linetype = 2,size=1.2)
p <- p + stat_summary(aes(x=mainAge,y=mv), fun=mean, colour="red", geom="point",group=1)
print(p)





##############
##  Fig 2A  ## # Proportion of Cardiomyocyte vs. hESC deregulated genes
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
#?

### Load datasets
df=read.table(paste0(datasets_dir,"CM_vs_hESC_DE_genes_Choy.txt"),header = T) # CM vs hESC Gene Exp

### Get deregulated genes (log2FC > 2 & padj < 0.1)
up=df[which(df$log2FoldChange>2 & df$padj<0.1),]
down=df[which(df$log2FoldChange<(-2) & df$padj<0.1),]

# dereg genes barplot - 16ages (supp?)
# up_vec=vector();down_vec=vector()
# for (i in 1:16){
#   up_vec=c(up_vec,nrow(up[which(up$Phylostrata==i),]))
#   down_vec=c(down_vec,nrow(down[which(down$Phylostrata==i),]))
# }
# tab=rbind(down_vec,up_vec);colnames(tab)=paste0("Class",1:16)

### Compute percentages of deregulated genes by main ages
up_vec_3ages=vector();down_vec_3ages=vector()
for (i in c("UC","EM","MM")){
  up_vec_3ages=c(up_vec_3ages,nrow(up[which(up$mainAge==i),]))
  down_vec_3ages=c(down_vec_3ages,nrow(down[which(down$mainAge==i),]))
}
tab_3ages=rbind(down_vec_3ages,up_vec_3ages);colnames(tab_3ages)=c("UC","EM","MM")
tab_3ages=rbind(tab_3ages,table(df$mainAge)[c(3,1,2)])
tab_3ages[1,]=tab_3ages[1,]/tab_3ages[3,]*100
tab_3ages[2,]=tab_3ages[2,]/tab_3ages[3,]*100
tab_3ages[3,]=100-(tab_3ages[1,]+tab_3ages[2,])

### Make Barplot
b=barplot(tab_3ages,col=c("red3","green3","grey"), ylab="",main="CM vs hESC deregulated genes \n Trigos Gene Ages",ylim=c(0,100),xaxt="n",axes=0,border = "white") #3mainAges
text(b, y= tab_3ages[1,]/2, labels = paste0(round(tab_3ages[1,]),"%"),cex=1.4)
text(b, y= tab_3ages[1,]+tab_3ages[2,]/2, labels = paste0(round(tab_3ages[2,]),"%"),cex=1.4)
text(b, y= tab_3ages[1,]+tab_3ages[2,]+tab_3ages[3,]/2, labels = paste0(round(tab_3ages[3,]),"%"),cex=1.4)
mtext("Genes (%)", side = 2, line = 2.5,cex=1.4)
axis(2, at=c(0,25,50,75,100),las=2,cex.axis=1.2,lwd=1.2,font=1)
text(cex=1.6, x=b+.4, y=-12, colnames(tab_3ages), xpd=T,pos=2)





##############
##  Fig 2B  ## # Number of Polycomb target genes in hESC, B-Cell and B-CLL
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75) #hg19
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
pc_hESC=read.table(paste0(datasets_dir,"H3K27me3_H9_hESC_WT_Yan_hg19_peaks.broadPeak"),header = F) # hESC H3K27me3 peaks
pc_bcel=read.table(paste0(datasets_dir,"H3K27me3_BCell_S0040EH1.ERX550664.20150527_broad_hg19.bed"),header = F) # B-Cell H3K27me3 peaks
pc_cll=read.table(paste0(datasets_dir,"H3K27me3_CLL_S00AYXH1.ERX604488.20150527_broad_hg19.bed"),header = F) # B-CLL H3K27me3 peaks

### Get gene TSS coordinates (from TSS-1bp to TSS+1bp)
genes=genes(EnsDb.Hsapiens.v75)
genes=genes[which(genes$gene_name %in% trigos$GeneID)]
seqlevels(genes)=paste0("chr",seqlevels(genes))
tss=genes
start=ifelse(strand(tss)=="+",start(tss)-1,end(tss)-1)
end=ifelse(strand(tss)=="+",start(tss)+1,end(tss)+1)
start(tss)=start;end(tss)=end

### Get genes which TSS overlaps with H3K27me3 peaks in hESC
pc_hESC=GRanges(pc_hESC$V1,IRanges(pc_hESC$V2,pc_hESC$V3))
ol=findOverlaps(tss,pc_hESC)
hESC_genes=genes[unique(queryHits(ol))]
hESC_genes=trigos[which(trigos$GeneID %in% hESC_genes$gene_name),]
hESC_genes=unique(hESC_genes[,c(1,4,5)])

### Get genes which TSS overlaps with H3K27me3 peaks in B-Cell
pc_bcel=GRanges(pc_bcel$V1,IRanges(pc_bcel$V2,pc_bcel$V3))
ol=findOverlaps(tss,pc_bcel)
bcel_genes=genes[unique(queryHits(ol))]
bcel_genes=trigos[which(trigos$GeneID %in% bcel_genes$gene_name),]
bcel_genes=unique(bcel_genes[,c(1,4,5)])

### Get genes which TSS overlaps with H3K27me3 peaks in B-CLL
pc_cll=GRanges(pc_cll$V1,IRanges(pc_cll$V2,pc_cll$V3))
ol=findOverlaps(tss,pc_cll)
cll_genes=genes[unique(queryHits(ol))]
cll_genes=trigos[which(trigos$GeneID %in% cll_genes$gene_name),]
cll_genes=unique(cll_genes[,c(1,4,5)])

### Make Barplot
color_3cells=c("#f69822","#598dc7","#b1332c")
m=c(nrow(hESC_genes),nrow(bcel_genes),nrow(cll_genes));names(m)=c("hESC","B-Cell","B-CLL")
b=barplot(m,col=c("#f69822","#598dc7","#b1332c"), ylab="",main="Polycomb target genes \n Trigos Gene Ages - TSS",ylim=c(0,4000),xaxt="n",axes=0,border = "white")
pos=ifelse(m>0,round(m,1),"")
text(b, m, labels=pos, cex= 1.6,pos=3)
mtext("Polycomb target genes", side = 2, line = 2,cex=1.4)
axis(2, at=c(0,4000),las=2,cex.axis=1.4,lwd=1.2,font=1)
text(cex=1.6, x=b+.3, y=-250, names(m), xpd=T, srt=45,pos=2)





##############
##  Fig 2C  ## # Percentages of Polycomb target genes by age in hESC, B-Cell and B-CLL
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75) #hg19
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
pc_hESC=read.table(paste0(datasets_dir,"H3K27me3_H9_hESC_WT_Yan_hg19_peaks.broadPeak"),header = F) # hESC H3K27me3 peaks
pc_bcel=read.table(paste0(datasets_dir,"H3K27me3_BCell_S0040EH1.ERX550664.20150527_broad_hg19.bed"),header = F) # B-Cell H3K27me3 peaks Blueprint
pc_cll=read.table(paste0(datasets_dir,"H3K27me3_CLL_S00AYXH1.ERX604488.20150527_broad_hg19.bed"),header = F) # B-CLL H3K27me3 peaks

### Get gene TSS coordinates (from TSS-1bp to TSS+1bp)
genes=genes(EnsDb.Hsapiens.v75)
genes=genes[which(genes$gene_name %in% trigos$GeneID)]
seqlevels(genes)=paste0("chr",seqlevels(genes))
tss=genes
start=ifelse(strand(tss)=="+",start(tss)-1,end(tss)-1)
end=ifelse(strand(tss)=="+",start(tss)+1,end(tss)+1)
start(tss)=start;end(tss)=end

### Get genes which TSS overlaps with H3K27me3 peaks in hESC
pc_hESC=GRanges(pc_hESC$V1,IRanges(pc_hESC$V2,pc_hESC$V3))
ol=findOverlaps(tss,pc_hESC)
hESC_genes=genes[unique(queryHits(ol))]
hESC_genes=trigos[which(trigos$GeneID %in% hESC_genes$gene_name),]
hESC_genes=unique(hESC_genes[,c(1,4,5)])

### Get genes which TSS overlaps with H3K27me3 peaks in B-Cell
pc_bcel=GRanges(pc_bcel$V1,IRanges(pc_bcel$V2,pc_bcel$V3))
ol=findOverlaps(tss,pc_bcel)
bcel_genes=genes[unique(queryHits(ol))]
bcel_genes=trigos[which(trigos$GeneID %in% bcel_genes$gene_name),]
bcel_genes=unique(bcel_genes[,c(1,4,5)])

### Get genes which TSS overlaps with H3K27me3 peaks in B-CLL
pc_cll=GRanges(pc_cll$V1,IRanges(pc_cll$V2,pc_cll$V3))
ol=findOverlaps(tss,pc_cll)
cll_genes=genes[unique(queryHits(ol))]
cll_genes=trigos[which(trigos$GeneID %in% cll_genes$gene_name),]
cll_genes=unique(cll_genes[,c(1,4,5)])

### Get numbers of polycomb target genes by age - 3 ages
m_hESC=vector();m_bcel=vector();m_cll=vector()
for(i in c("UC","EM","MM")){
  m_hESC=c(m_hESC,nrow(hESC_genes[which(hESC_genes$mainAge==i),]))
  m_bcel=c(m_bcel,nrow(bcel_genes[which(bcel_genes$mainAge==i),]))
  m_cll=c(m_cll,nrow(cll_genes[which(cll_genes$mainAge==i),]))
}
names(m_hESC)=c("UC","EM","MM");names(m_bcel)=c("UC","EM","MM");names(m_cll)=c("UC","EM","MM")

### Get percentages of polycomb target genes by age - 3 ages
m_hESC=m_hESC/table(trigos$mainAge)[c(3,1,2)]*100;m_hESC=rbind(m_hESC,c(50-m_hESC))
m_bcel=m_bcel/table(trigos$mainAge)[c(3,1,2)]*100;m_bcel=rbind(m_bcel,c(50-m_bcel))
m_cll=m_cll/table(trigos$mainAge)[c(3,1,2)]*100;m_cll=rbind(m_cll,c(50-m_cll))

### Make Barplots
for (cell in c("hESC","B-Cell","B-CLL")){
  if(cell=="hESC"){tab_3ages=m_hESC;col="#f69822"}
  if(cell=="B-Cell"){tab_3ages=m_bcel;col="#598dc7"}
  if(cell=="B-CLL"){tab_3ages=m_cll;col="#b1332c"}
b=barplot(tab_3ages,col=c(col,"grey"), ylab="",main=paste0(cell," - Polycomb target genes - TSS ol"),ylim=c(0,50),xaxt="n",axes=0,border = "white")
text(b, y= tab_3ages[1,], labels = paste0(round(tab_3ages[1,]),"%"),cex=1.6,pos=3)
mtext("Polycomb target genes (%)", side = 2, line = 2.5,cex=1.4)
axis(2, at=c(0,25,50),las=2,cex.axis=1.4,lwd=1.2,font=1)
text(cex=2, x=b, y=-7, colnames(tab_3ages), xpd=T)
}





##############
##  Fig 2E  ## # PolII Pausing Index by age in hESC, B-Cell and B-CLL
##############
### N.B: PolII Pausing Indexes computed with script : "Paused_PolII_Index_computing.R"


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(ggpubr)
#?

### Load datasets
trigos=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
load(paste0(datasets_dir,"hESC_Paused_PolII_Index_table.Rdata")) #hESC_index
load(paste0(datasets_dir,"B-Cell_Paused_PolII_Index_table.Rdata")) #bcel_index
load(paste0(datasets_dir,"Lymphoma_Paused_PolII_Index_table.Rdata")) #lymphoma_index

### Remove NA and Inf values
#hESC
hESC_index=hESC_index[which(!(is.na(hESC_index$index))),]
hESC_index=hESC_index[which(!(is.infinite(hESC_index$index))),]
hESC_index=unique(hESC_index)
hESC_index$Cell="hESC"

#B-Cell
bcel_index=bcel_index[which(!(is.na(bcel_index$index))),]
bcel_index=bcel_index[which(!(is.infinite(bcel_index$index))),]
bcel_index=unique(bcel_index)
bcel_index$Cell="B-Cell"

#B-Lymphoma
lymphoma_index=lymphoma_index[which(!(is.na(lymphoma_index$index))),]
lymphoma_index=lymphoma_index[which(!(is.infinite(lymphoma_index$index))),]
lymphoma_index=unique(lymphoma_index)
lymphoma_index$Cell="Lymphoma"

### Remove outliers
#hESC
Q1 <- quantile(hESC_index$index, .25);Q3 <- quantile(hESC_index$index, .75);IQR <- IQR(hESC_index$index)
no_outliers <- subset(hESC_index, hESC_index$index > (Q1 - 1.5*IQR) & hESC_index$index < (Q3 + 1.5*IQR))
hESC_index=no_outliers
#B-Cell
Q1 <- quantile(bcel_index$index, .25);Q3 <- quantile(bcel_index$index, .75);IQR <- IQR(bcel_index$index)
no_outliers <- subset(bcel_index, bcel_index$index > (Q1 - 1.5*IQR) & bcel_index$index < (Q3 + 1.5*IQR))
bcel_index=no_outliers
#Lymphoma
Q1 <- quantile(lymphoma_index$index, .25);Q3 <- quantile(lymphoma_index$index, .75);IQR <- IQR(lymphoma_index$index)
no_outliers <- subset(lymphoma_index, lymphoma_index$index > (Q1 - 1.5*IQR) & lymphoma_index$index < (Q3 + 1.5*IQR))
lymphoma_index=no_outliers

colnames(hESC_index)[5]="index_hESC";colnames(CM_index)[5]="index_CM"
colnames(bcel_index)[5]="index_bcel";colnames(lymphoma_index)[5]="index_lymphoma"
df=merge(hESC_index[,c(1,2,5)],bcel_index[,c(1,5)],by="Gene")
df=merge(df,lymphoma_index[,c(1,5)],by="Gene")

### Remove short genes (length < 1kb)
genes=as.data.frame(genes(EnsDb.Hsapiens.v75))
genes$width=genes$end-genes$start
remove=genes[which(genes$width<1000),7]
df=df[which(!(df$Gene %in% remove)),]

### Remove unexpressed genes
keep <- rowSums(df[,3:5]) >= 0.2
df <- df[keep,]

### Remove multiple annotated genes
remove=table(df$Gene)
remove=names(remove[which(remove>1)])
df=df[which(!(df$Gene %in% remove)),]

### Normalization between 0 and 1
df$index_hESC_norm=(df$index_hESC-min(df$index_hESC))/(max(df$index_hESC)-min(df$index_hESC))
df$index_bcel_norm=(df$index_bcel-min(df$index_bcel))/(max(df$index_bcel)-min(df$index_bcel))
df$index_lymphoma_norm=(df$index_lymphoma-min(df$index_lymphoma))/(max(df$index_lymphoma)-min(df$index_lymphoma))

### Group by main age
df$mainAge=NA
df$mainAge=ifelse(df$Phylostrata<4,"UC",NA)
df$mainAge=ifelse(df$Phylostrata>3 & df$Phylostrata<11,"EM",df$mainAge)
df$mainAge=ifelse(df$Phylostrata>10,"MM",df$mainAge)

### Make Boxplot
color_3cells=c("#f69822","#598dc7","#b1332c")
my_comparisons <- list(c("hESC","BCell") ,c("BCell","Lymphom"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
df2=data.frame(Gene=rep(df$Gene,3),Phylostrata=rep(df$Phylostrata,3),mainAge=rep(df$mainAge,3),Index_norm=c(df$index_hESC_norm,df$index_bcel_norm,df$index_lymphoma_norm),Cell=c(rep("hESC",nrow(df)),rep("BCell",nrow(df)),rep("Lymphom",nrow(df))))
df2$mainAge=factor(df2$mainAge,levels = c("UC","EM","MM"))

### Option 1 : Sort by main age
# bxp <- ggboxplot(df2, x = "Cell", y = "Index_norm",
#                  facet.by = "mainAge",fill="Cell",palette=color_3cells,paired=T) +ylim(0,1.2) + ylab("Pausing Index") + xlab("")+
#   stat_compare_means(comparisons = my_comparisons,method="wilcox.test",symnum.args=symnum.args)
# bxp

### Option 2 : Sort by cell type  
df2$Cell=factor(df2$Cell,levels = c("hESC","BCell","Lymphom"))
my_comparisons <- list(c("EM","MM"),c("UC","EM"))
bxp <- ggboxplot(df2, x = "mainAge", y = "Index_norm",
                 facet.by = "Cell",fill="Cell",palette=color_3cells,paired=T) +ylim(0,1.2) + ylab("Pausing Index") + xlab("")+
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",symnum.args=symnum.args) 
bxp





##############
##  Fig 3A  ## # Percentages of LAD genes by age
##############


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
#?

### Load datasets
ages=read.table(paste0(datasets_dir,"Gene_ages_Trigos.txt"),header = T) # Gene Ages
lad=read.table(paste0(datasets_dir,"LADatlas_100kbBin_LP160718_hg19.bed"),header = F) # LADs dataset

### Get age-annotated gene promoters coordinates (from TSS-30bp to TSS+300bp)
prom=genes(EnsDb.Hsapiens.v75)
seqlevels(prom)=paste0("chr",seqlevels(prom))
start=ifelse(strand(prom)=="+",start(prom)-30,end(prom)-300)
end=ifelse(strand(prom)=="+",start(prom)+300,end(prom)+30)
start(prom)=start;end(prom)=end

### Get LAD overlapping genes
lad=lad[which(lad$V6=="fLAD" | lad$V6=="cLAD"),]
lad=GRanges(lad$V1,IRanges(lad$V2,lad$V3))
ol=findOverlaps(prom,lad)
lad_genes=prom[queryHits(ol)]$gene_id
lad_genes=ages[which(ages$gene_id %in% lad_genes),]

### Get percentages
m=table(lad_genes$mainAge)[c(3,1,2)]
m=rbind(m,table(ages$mainAge)[c(3,1,2)])
m[1,]=m[1,]/m[2,]*100
m[2,]=20-(m[1,])

### Make Barplot
col=c("#ffb09c","grey")
b=barplot(m,col=col, ylab="",main="LADs Genes",ylim=c(0,20),xaxt="n",axes=0,border = "white") #3mainAges
text(x=b, y= m[1,], pos=3,labels = paste0(round(m[1,]),"%"),cex=1.4)
mtext("Genes (%)", side = 2, line = 2.5,cex=1.4)
axis(2, at=c(0,20),las=2,cex.axis=1.2,lwd=1.2,font=1)
text(cex=1.6, x=b+.3, y=-2, colnames(m), xpd=T,pos=2)





##############
##  Fig 3C  ## # Global ChAs of age-specific gene promoters in hESC, B-Cell and B-CLL
##############

datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
library(chaser)
#?

### Load datasets
hESC_PP=read.table(paste0(datasets_dir,"PCHiC_hESC_PP_Narrow_hg19.txt"),header = T) #hESC 'Narrow' PP
bcel_PP=read.table(paste0(datasets_dir,"PCHiC_BCell_PP_Narrow_hg19.txt"),header = T) #B-Cell 'Narrow' PP
cll_PP=read.table(paste0(datasets_dir,"PCHiC_BCLL_PP_Narrow_hg19.txt"),header = T) #B-CLL 'Narrow' PP 
feat=read.table(paste0(datasets_dir,"CHASER_Features_table_mainAges_ALL_networks.txt"),header = T)

### Get PP Networks
hESC_PP_net=make_chromnet(hESC_PP[,c(1:6)])
bcel_PP_net=make_chromnet(bcel_PP[,c(1:6)])
cll_PP_net=make_chromnet(cll_PP[c(1:6)])

net=hESC_PP_net
net=load_features(net,feat,type = "features_table",missingv = NA,featnames = "mainAge")

#turn into character 
for(i in 1:ncol(net$features)){
  net$features[,i]=as.character(net$features[,i])
}

#compute ChAs 
chas_net <- chas(net,"categorical")

#compute Randomizations
random_hic_net=chaser::randomize(net = net,nrandom = 100,dist.match = T)

nb_random=length(random_hic_net)

rm(chas_random_hic_net)
chas_random_hic_net=vector()
for (i in 1:length(random_hic_net)){
  chas_random_hic_net=c(chas_random_hic_net,chas(random_hic_net[[i]],"categorical")) #for categorical ChAs
}
value=chas_random_hic_net
sample=names(chas_random_hic_net)
rm(df_random_hic_net)
df_random_hic_net=data.frame()
df_random_hic_net=matrix(c(value,sample),length(chas_random_hic_net),ncol=2)
colnames(df_random_hic_net)=c("Value","Sample")
df_random_hic_net=as.data.frame(df_random_hic_net)
df_random_hic_net$Value=as.numeric(as.character(df_random_hic_net$Value))
value=chas_net
sample=names(chas_net)

rm(df_chas_hic_net)
df_chas_hic_net=data.frame()
df_chas_hic_net=matrix(c(value,sample),length(chas_net),ncol=2)
colnames(df_chas_hic_net)=c("Value","Sample")
df_chas_hic_net=as.data.frame(df_chas_hic_net)
df_chas_hic_net$Value=as.numeric(as.character(df_chas_hic_net$Value))
df_hic_net=rbind(df_random_hic_net,df_chas_hic_net)

####compute the z-score
nb_features=length(unique(df_hic_net$Sample))
test_random=df_hic_net[1:(nb_features*nb_random),]
test_chas=df_hic_net[-(1:(nb_features*nb_random)),]
mean_random <- aggregate(. ~ Sample, data = test_random, mean)
sd_random <- aggregate(. ~ Sample, data = test_random, sd)

rm(i,zscore)
zscore=vector()
for(i in unique(df_hic_net$Sample)){
  zscore=c(zscore,(test_chas[which(test_chas$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
}
names(zscore)=unique(df_hic_net$Sample)

### Make Z-Scores Barplot
color_3ages=c("#f7766d","#02bb38","#609bff")
zscore=zscore[c(3,1,2)]

## get published results
#zscore=readRDS(paste0(datasets_dir,"hESC_ages_zscores.rds"))  # hESC
#zscore=readRDS(paste0(datasets_dir,"BCell_ages_zscores.rds")) # B-Cell
#zscore=readRDS(paste0(datasets_dir,"CLL_ages_zscores.rds"))   # CLL

b=barplot(zscore,col=color_3ages, ylab="",
          main=paste0("ChAs Z-Scores \n","Network: true ann CLL PP"),
          ylim=c(0,10),xaxt="n",axes=0,border = "white") 
pos=ifelse(zscore>0,round(zscore,1),"")
neg=ifelse(zscore<0,round(zscore,1),"")
text(b, zscore, labels=pos, cex= 1.6,pos=3)
text(b, zscore, labels=neg, cex= 1.6,pos=1)
mtext("ChAs Z-Scores", side = 2, line = 2,cex=1.6)
axis(2, at=c(0,12),las=2,cex.axis=1.2,lwd=1.2,font=1)
text(cex=1.6, x=b+.4, y=-2, names(zscore), xpd=T,pos=2)





##############
##  Fig 3F  ## # ChAs of PolII, Polycomb, and LAD genes vs others in hESC, B-Cell, and B-CLL
##############
### N.B: DNA frag annotation file completed with script : "New_true_annotation_geneAges.R"


datasets_dir="/Users/fla/Desktop/GeneAges_Proper_Datasets/"


### Load libraries
library(EnsDb.Hsapiens.v75)
library(chaser)
#?

### Load datasets
hESC_PP=read.table(paste0(datasets_dir,"PCHiC_hESC_PP_Narrow_hg19.txt"),header = T) #hESC 'Narrow' PP
bcel_PP=read.table(paste0(datasets_dir,"PCHiC_BCell_PP_Narrow_hg19.txt"),header = T) #B-Cell 'Narrow' PP
cll_PP=read.table(paste0(datasets_dir,"PCHiC_BCLL_PP_Narrow_hg19.txt"),header = T) #B-CLL 'Narrow' PP 
ann=read.table(paste0(datasets_dir,"HindIII_annotation_ens37_completed.txt"),header = T) # DNA fragment annotation

### Get PP Networks
hESC_PP_net=make_chromnet(hESC_PP[,c(1:6)])
bcel_PP_net=make_chromnet(bcel_PP[,c(1:6)])
cll_PP_net=make_chromnet(cll_PP[c(1:6)])

net=hESC_PP_net

#gene ages + LAD genes
net=load_features(net,ann[which(ann$LAD_gene==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_LAD")
net=load_features(net,ann[which(ann$LAD_gene==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_LAD")
net=load_features(net,ann[which(ann$LAD_gene==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_LAD")
net=load_features(net,ann[which(ann$LAD_gene==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_LAD")
net=load_features(net,ann[which(ann$LAD_gene==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_LAD")
net=load_features(net,ann[which(ann$LAD_gene==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_LAD")
colnames(net$features)=c("UC_LAD","UC_no_LAD","EM_LAD","EM_no_LAD","MM_LAD","MM_no_LAD")

#gene ages + Pc TSS genes
#hESC
net=load_features(net,ann[which(ann$hESC_Pc_TSS==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_Pc")
net=load_features(net,ann[which(ann$hESC_Pc_TSS==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_Pc")
net=load_features(net,ann[which(ann$hESC_Pc_TSS==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_Pc")
net=load_features(net,ann[which(ann$hESC_Pc_TSS==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_Pc")
net=load_features(net,ann[which(ann$hESC_Pc_TSS==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_Pc")
net=load_features(net,ann[which(ann$hESC_Pc_TSS==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_Pc")
colnames(net$features)=c("UC_Pc","UC_no_Pc","EM_Pc","EM_no_Pc","MM_Pc","MM_no_Pc")
#B-Cell
net=load_features(net,ann[which(ann$bcel_Pc_TSS==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_Pc")
net=load_features(net,ann[which(ann$bcel_Pc_TSS==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_Pc")
net=load_features(net,ann[which(ann$bcel_Pc_TSS==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_Pc")
net=load_features(net,ann[which(ann$bcel_Pc_TSS==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_Pc")
net=load_features(net,ann[which(ann$bcel_Pc_TSS==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_Pc")
net=load_features(net,ann[which(ann$bcel_Pc_TSS==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_Pc")
colnames(net$features)=c("UC_Pc","UC_no_Pc","EM_Pc","EM_no_Pc","MM_Pc","MM_no_Pc")
#CLL
net=load_features(net,ann[which(ann$cll_Pc_TSS==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_Pc")
net=load_features(net,ann[which(ann$cll_Pc_TSS==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_Pc")
net=load_features(net,ann[which(ann$cll_Pc_TSS==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_Pc")
net=load_features(net,ann[which(ann$cll_Pc_TSS==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_Pc")
net=load_features(net,ann[which(ann$cll_Pc_TSS==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_Pc")
net=load_features(net,ann[which(ann$cll_Pc_TSS==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_Pc")
colnames(net$features)=c("UC_Pc","UC_no_Pc","EM_Pc","EM_no_Pc","MM_Pc","MM_no_Pc")

#gene ages + PolII TSS genes
#hESC
net=load_features(net,ann[which(ann$hESC_Pol_TSS==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_pol")
net=load_features(net,ann[which(ann$hESC_Pol_TSS==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_pol")
net=load_features(net,ann[which(ann$hESC_Pol_TSS==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_pol")
net=load_features(net,ann[which(ann$hESC_Pol_TSS==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_pol")
net=load_features(net,ann[which(ann$hESC_Pol_TSS==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_pol")
net=load_features(net,ann[which(ann$hESC_Pol_TSS==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_pol")
colnames(net$features)=c("UC_pol","UC_no_pol","EM_pol","EM_no_pol","MM_pol","MM_no_pol")
#B-Cell
net=load_features(net,ann[which(ann$bcel_Pol_TSS==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_pol")
net=load_features(net,ann[which(ann$bcel_Pol_TSS==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_pol")
net=load_features(net,ann[which(ann$bcel_Pol_TSS==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_pol")
net=load_features(net,ann[which(ann$bcel_Pol_TSS==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_pol")
net=load_features(net,ann[which(ann$bcel_Pol_TSS==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_pol")
net=load_features(net,ann[which(ann$bcel_Pol_TSS==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_pol")
colnames(net$features)=c("UC_pol","UC_no_pol","EM_pol","EM_no_pol","MM_pol","MM_no_pol")
#CLL
net=load_features(net,ann[which(ann$cll_Pol_TSS==2 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_pol")
net=load_features(net,ann[which(ann$cll_Pol_TSS==1 & grepl("UC",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "UC_no_pol")
net=load_features(net,ann[which(ann$cll_Pol_TSS==2 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_pol")
net=load_features(net,ann[which(ann$cll_Pol_TSS==1 & grepl("EM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "EM_no_pol")
net=load_features(net,ann[which(ann$cll_Pol_TSS==2 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_pol")
net=load_features(net,ann[which(ann$cll_Pol_TSS==1 & grepl("MM",ann$main_Ages)),c(1:3,12)],type = "features_table",missingv = NA,featnames = "MM_no_pol")
colnames(net$features)=c("UC_pol","UC_no_pol","EM_pol","EM_no_pol","MM_pol","MM_no_pol")

#Process matrix for categorical ChAs
for(i in 1:ncol(net$features)){
  net$features[,i]=ifelse(net$features[,i]>0,2,1)
  net$features[,i]=ifelse(is.na(net$features[,i]),1,net$features[,i])
  net$features[,i]=as.character(net$features[,i])
}

chas_net=chas(net,"categorical")

random_hic_net=chaser::randomize(net = net,nrandom = 100,dist.match = T)

nb_random=length(random_hic_net)

rm(chas_random_hic_net)
chas_random_hic_net=vector()
for (i in 1:length(random_hic_net)){
  chas_random_hic_net=c(chas_random_hic_net,chas(random_hic_net[[i]],"categorical")) #for categorical ChAs
}
value=chas_random_hic_net
sample=names(chas_random_hic_net)
rm(df_random_hic_net)
df_random_hic_net=data.frame()
df_random_hic_net=matrix(c(value,sample),length(chas_random_hic_net),ncol=2)
colnames(df_random_hic_net)=c("Value","Sample")
df_random_hic_net=as.data.frame(df_random_hic_net)
df_random_hic_net$Value=as.numeric(as.character(df_random_hic_net$Value))
value=chas_net
sample=names(chas_net)

rm(df_chas_hic_net)
df_chas_hic_net=data.frame()
df_chas_hic_net=matrix(c(value,sample),length(chas_net),ncol=2)
colnames(df_chas_hic_net)=c("Value","Sample")
df_chas_hic_net=as.data.frame(df_chas_hic_net)
df_chas_hic_net$Value=as.numeric(as.character(df_chas_hic_net$Value))
df_hic_net=rbind(df_random_hic_net,df_chas_hic_net)

####compute the z-score
nb_features=length(unique(df_hic_net$Sample))
test_random=df_hic_net[1:(nb_features*nb_random),]
test_chas=df_hic_net[-(1:(nb_features*nb_random)),]
mean_random <- aggregate(. ~ Sample, data = test_random, mean)
sd_random <- aggregate(. ~ Sample, data = test_random, sd)

rm(i,zscore)
zscore=vector()
for(i in unique(df_hic_net$Sample)){
  zscore=c(zscore,(test_chas[which(test_chas$Sample==i),1]-mean_random[which(mean_random$Sample==i),2])/sd_random[which(sd_random$Sample==i),2])
}
names(zscore)=unique(df_hic_net$Sample)

delta_col=c(rep("#f7766d",2),rep("#02bb38",2),rep("#609bff",2))

b=barplot(zscore,col=delta_col, ylab="",
          main="PP Narrow hESC Network\nNew True ann",
          ylim=c(0,15),xaxt="n",axes=0,border = "white") 
pos=ifelse(zscore>0,round(zscore,2),"")
neg=ifelse(zscore<0,round(zscore,2),"")
text(b, zscore, labels=pos, cex= 1.4,pos=3)
text(b, zscore, labels=neg, cex= 1.4,pos=1)
mtext("ChAs Z-Scores", side = 2, line = 2,cex=1.6)
axis(2, at=c(0,15),las=2,cex.axis=1.4,lwd=1.2,font=1)
text(cex=1.2, x=b+.4, y=-1, names(zscore), xpd=T,pos=2,srt=45)

#get published results
# zscore=readRDS(paste0(datasets_dir,"hESC_PolII_zscores.rds")) # hESC PolII
# zscore=readRDS(paste0(datasets_dir,"bcel_PolII_zscores.rds")) # B-Cell PolII
# zscore=readRDS(paste0(datasets_dir,"cll_PolII_zscores.rds")) # CLL PolII
# zscore=readRDS(paste0(datasets_dir,"hESC_Pc_zscores.rds")) # hESC Pc
# zscore=readRDS(paste0(datasets_dir,"bcel_Pc_zscores.rds")) # B-Cell Pc
# zscore=readRDS(paste0(datasets_dir,"cll_Pc_zscores.rds")) # CLL Pc
# zscore=readRDS(paste0(datasets_dir,"hESC_LAD_zscores.rds")) # hESC LAD
# zscore=readRDS(paste0(datasets_dir,"bcel_LAD_zscores.rds")) # B-Cell LAD
# zscore=readRDS(paste0(datasets_dir,"cll_LAD_zscores.rds")) # CLL LAD

#Delta diff
diff=c(zscore[1]-zscore[2],zscore[3]-zscore[4],zscore[5]-zscore[6])
names(diff)=c("UC","EM","MM")
col=c("#f7766d","#02bb38","#609bff")

b=barplot(diff,col=col, ylab="",
          main="Delta ChAs - LAD genes\nPP Narrow CLL Network - New True ann",
          ylim=c(-15.8,15.8),xaxt="n",axes=0,border = "white") 
pos=ifelse(diff>0,round(diff,1),"")
neg=ifelse(diff<0,round(diff,1),"")
text(b, diff, labels=pos, cex= 1.6,pos=3)
text(b, diff, labels=neg, cex= 1.6,pos=1)
mtext("Δ ChAs", side = 2, line = 2,cex=1.8) 
axis(2, at=c(-15,0,15),las=2,cex.axis=1.4,lwd=1.2,font=1)
text(cex=1.6, x=b+.4, y=-20, names(diff), xpd=T,pos=2)




