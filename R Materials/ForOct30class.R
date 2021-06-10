install.packages("ggplot2")
install.packages("ploty")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
				
library(edgeR)
library(ggplot2)
library(plotly)

#getting CCNA numbers and descriptors
setwd("R Materials/")
CCNA_description <- read.delim("CCNA_description.txt")

#getting the kegg terms 
GK <- getGeneKEGGLinks(species.KEGG ="ccs") ##other species "dre" --> Zebra fish
GPN <- getKEGGPathwayNames(species.KEGG="ccs")
#newgk <- GK
#names(newgk) <- c("genename","PathwayID")


x <- read.table("Compiled_norRNA.txt",row.names="gene",header=T)
group <- factor(c(1,1,1,2,2,2,3,3,3,4,4,4))
y <- DGEList(counts=x,group=group)

keep <- rowSums(cpm(y)>30) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)

DLvWTqlf <- glmQLFTest(fit,coef=2)
E4vWTqlf <- glmQLFTest(fit,coef=3)
XGvWTqlf <- glmQLFTest(fit,coef=4)
XGvDLqlf <- glmQLFTest(fit,contrast=c(0,-1,0,1))


#showing kegg terms
#important to change qlf file in these two lines and path name

qlffile <- DLvWTqlf

one <- topTags(qlffile,n=3864)
keg <- kegga(qlffile,species.KEGG="ccs",FDR=0.05)
topKEGG(keg,sort="up")
topKEGG(keg,sort="down")

#important to change KEGG path here to what you care about
kegpath <- "path:ccs03010"

#make a new data frame and mark those with FDR < 0.05
df <- one$table
df$signifcance <- ifelse(df$FDR<0.05,"FDR < 0.05","FDR > 0.05")
df$genename <- row.names(one$table)
newish <- merge(df,CCNA_description,by="genename")
new <- newish[order(newish$FDR),]
kegdescriptor <- GPN$Description[GPN$PathwayID==kegpath]

# extract the gene names from a certain KEGG path 
i <- row.names(one$table) %in% GK$GeneID[GK$PathwayID==kegpath]
newvec1 <- row.names(one$table[i,])

dfpath1 <- subset(new, (new$genename %in% newvec1))
dfpath1 <- subset(dfpath1,dfpath1$FDR<0.1)

#using ggplot2 to make the volcano plot for the data
volc <- ggplot(new,aes(x=logFC,y=-log10(PValue),text=paste(genename,description,sep="\n"),color=signifcance))+geom_point(alpha=0.6)+scale_color_manual(values=c("red","black"))
volc
volc + geom_point(data=subset(dfpath1),aes(x=logFC,y=-log10(PValue)),color="blue",size=4) + ggtitle(paste("Deletion-Lon v WT",kegdescriptor,sep="\n")) + theme_classic()

#uses gplotly to make an interactive plot
ggplotly()
