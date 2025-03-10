#####step1.bulk#######
library(IOBR)
library(easybio)
library(limma)
library(export)
setwd('/Users/drning/Desktop/ICH/1.bulk分析')
#cohort1
count=read.csv('exp_cohort1.csv')
count<-remove_duplicate_genes(eset = count, column_of_symbol = "GeneName", method = "mean")
tpm1=easy_count2data(data = count1, 
                     to = "tpm",  
                     idType = "Symbol", 
                     org = "hsa", 
                     needlog = T
)
tpm1=easy_filtergene(data = tpm1,
                     genecode = "protein_coding", 
                     org = "hs",
                     hgVersion = "v38", 
                     filtersd = 0, 
                     filterrowsum = 0 
)
#cohort2
count=read.csv('exp_cohort2.csv')
count<-remove_duplicate_genes(eset = count, column_of_symbol = "GeneName", method = "mean")
tpm2=easy_count2data(data = count, 
                     to = "tpm",  
                     idType = "Symbol", 
                     org = "hsa", 
                     needlog = T 
)
tpm2=easy_filtergene(data = tpm2,
                     genecode = "protein_coding", 
                     org = "hs", 
                     hgVersion = "v38", 
                     filtersd = 0, 
                     filterrowsum = 0 
)
gene=intersect(rownames(tpm1),rownames(tpm2))
tpm1=tpm1[gene,]
tpm2=tpm2[gene,]

boxplot(tpm1,las=2,cex.axis=0.6)
graph2pdf(file='cohort1boxplot.pdf',height=3,width=8)
boxplot(tpm2,las=2,cex.axis=0.6)
graph2pdf(file='cohort2boxplot.pdf',height=4,width=8)

tpm=bind_cols(tpm1,tpm2)
cli=read.table('/Users/drning/Desktop/ICH/1.bulk分析/group.txt', header=T, sep="\t", check.names=F, row.names=1)  
cli=cli[colnames(tpm),]
boxplot(tpm,las=2,cex.axis=0.6)
graph2pdf(file='Normbatch.boxplot.pdf',height=4,width=9)
PCAres=easy_PCAplot(expr = tpm,
                    scale = T,addcir = T,
                    Group =c(rep('Cohort1',92),rep('Cohort2',38)),
                    levels = c("Cohort1", "Cohort2"),
                    cols = c("#E64B35FF","#1F78B4" ))
graph2pdf(file='Normbatch.PCA.pdf',height=4,width=5)

library(sva)
tpm=ComBat(tpm, cli$cohort, par.prior=TRUE)
boxplot(tpm,las=2,cex.axis=0.6)
graph2pdf(file='rmbatch.boxplot.pdf',height=4,width=9)
PCAres=easy_PCAplot(expr = tpm,
                    scale = T,addcir = T,
                    Group = c(rep('Cohort1',92),rep('Cohort2',38)),
                    levels =  c("Cohort1", "Cohort2"),
                    cols = c("#E64B35FF","#1F78B4" ))
graph2pdf(file='rmbatch.PCA.pdf',height=4,width=5)

cli=arrange(cli,group)
tpm=tpm[,rownames(cli)]
PCAres=easy_PCAplot(expr = tpm,
                    scale = T,
                    Group = cli$group,
                    levels = c("HBP", "ICH"),
                    cols = c("#A6CEE3" ,"#1F78B4"))
graph2pdf(file='PCA1.pdf',height=4,width=5)
PCAres=easy_PCAplot(expr = tpm,
                    scale = T,
                    Group = cli$cohort,
                    levels = c("cohort1", "cohort2"),
                    cols = c("#3BBCA8", "#FDBF6F"))
graph2pdf(file='PCA2.pdf',height=4,width=5)

DEG=easy_DEGs(data=tpm, 
              countdata=F,
              groups=cli$group, 
              treatname='ICH', 
              controlname='HBP',
              log2FC=0.5, 
              Pmeth="FDR",
              Pcutoff=0.05, 
              outfile="DEGout.txt" 
)
easy_volcano1(DEG = DEG, 
             logFC.Ncol = 2,
             Select.P = "FDR",
             P.Ncol = 6, 
             DEG.type.Ncol = 8,
             cutoff.P = 0.05,
             cutoff.logFC = 0.5,
             colors = c( "#A6CEE3"  ,"#E3E3E3","#1F78B4"))
graph2pdf(file='Volcan.pdf',height=4,width=4)
save(tpm,file = 'tpm.Rdata')
#heatmap
library(ComplexHeatmap)
library(viridis)
library(circlize)
DEG=read.table('DEGout.txt', header=T, sep="\t", check.names=F)  
diffSig <- DEG
rt <- tpm
updiffSig <- filter(diffSig,Type %in% c('Up')) 
upexp=rt[updiffSig$Gene,]
DowndiffSig <- filter(diffSig,Type %in% c('Down')) 
downexp=rt[DowndiffSig$Gene,]
hmExp=bind_rows(downexp,upexp)
head(hmExp)
hmExp=t(scale(t(hmExp)))
cli=read.table('group.txt', header=T, sep="\t", check.names=F, row.names=1)  
my=cli[colnames(hmExp),]
group<-  c("#A6CEE3" ,"#1F78B4")
my$group<- factor(my$group,levels = c('HBP','ICH'))
names(group) <- levels(my$group)
cohort<-  c("#3BBCA8", "#FDBF6F")
my$cohort<- factor(my$cohort,levels = c('cohort1','cohort2'))
names(cohort) <- levels(my$cohort)
Top = HeatmapAnnotation(group=my$group,
                        cohort=my$cohort,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,
                                                     title_gp = gpar(fontsize = 10,fontface = "bold"),
                                                     ncol=1),
                        border = T,
                        col=list(group=group ,
                                 cohort=cohort),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))
label_gene <- c('MMP8','ARG1','PLAU','RAB27A',#neutrophil granulocyte degranulation
                     'FCGR1A','HMGB2','S100A12','TLR5','IL18R1','IL18',#inflammatory response
                     'CD247','CD8A','CD8B','CCL5','ZAP70','NKG7',#lymphocyte activation
                     'PRF1','LCK','LAT',#PID CD8 TCR PATHWAY
                     'SMPD3','PDGFB','PDGFRB'#nuclear division
)  
Heatmap(hmExp, name = 'Z-score',
        top_annotation = Top,
        cluster_rows = F, show_row_names = F,
        col = colorRamp2(c(-2, 0, 2), c("#3E94B5", "white", "#ED6355")),
        color_space = "RGB",
        cluster_columns = FALSE, border = T,
        row_order = NULL,
        row_names_side = 'left',
        column_order = NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 8),
        column_split = c(rep(1, 64), rep(2, 66)),
        row_split = c(rep(1, 825), rep(2, 645)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 10), border = T,
                                    title_gp = gpar(fontsize = 10, fontface = "bold")),
        column_gap = unit(2, 'mm'),
        row_names_max_width = unit(22, 'cm'),
        left_annotation = rowAnnotation(
          mark_gene = anno_mark(
            at = which(rownames(hmExp) %in% label_gene), 
            labels = label_gene,  
            labels_gp = gpar(fontsize = 4)  
          )
        ),use_raster = T
)
graph2pdf(file='heatmapDEG.pdf',height=3,width=4)
######step2.scRNA########
library(Seurat)
library(DT)
library(COSG)
library(dplyr)
library(scRNAtoolVis)
library(data.table)
library(readr)
library(dplyr)
library(plyr)
library(ggpubr)
library(ggplot2)
library(patchwork)
mycolor <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
             "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#4DBBD5FF", 
             "#E64B35FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", 
             "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#BC3C29FF", "#0072B5FF", 
             "#E18727FF", "#D24B27", "#D51F26", "#272E6A", "#208A42", "#89288F", 
             "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", "#90D5E4", 
             "#89C75F", "#F37B7D", "#9983BD", "#3BBCA8", "#6E4B9E", "#0C727C", 
             "#7E1416", "#D8A767", "#3D3D3D")
setwd('/Users/drning/Desktop/ICH/2.scRNA注释')
load('mousesce.Rdata')
Idents(sce) <- 'celltype'
ann=c('Neutrophil', 'Monocyte', 'T', 'B', 'Megakaryocytes', 'Erythroblasts', 'NK')
sce$celltype <- factor(sce$celltype,levels = ann)
gene=c('Clec4d', 'Csf3r', 'Cxcr2',
       'Ccr2',  'Cybb', 'Klf4', 
       'Cd3g','Cd3d','Cd3e',
       'Cd79a', 'Cd79b', 'Cd19', 
       'Pf4','Cd151', 'Gp1ba',
       'Ctsw', 'Gzma',  'Klrd1')
jjDotPlot(object = sce,assay='RNA',slot='data',
          gene = gene,  id = 'celltype', 
          split.by.aesGroup =T, 
          xtree = F,ytree = F,
          rescale = T,rescale.min = -2,rescale.max =2,   midpoint = 0, dot.col = c('grey','white','#E18727FF'),
          point.geom = T,point.shape = 21,tile.geom = F)
graph2pdf(file='mousedotplotann.pdf',width=8,height=6)
SCP::CellDimPlot(
  srt = sce, group.by = "celltype",
  reduction = "UMAP", bg_color = "white",
  pt.size = 2,pt.alpha = 1,raster = T,
  palcolor=mycolor,
  split.by = 'group',
  label = F,label.size =4,label_insitu = TRUE,label_repel = F,
  add_density = F, 
  theme_use = "theme_blank"
)
graph2pdf(file='mousedim.pdf',width=8,height=6)


meta.tb = sce@meta.data
dat.plot = data.table(meta.tb)
dat.plot$Celltype = dat.plot$celltype
dat.plot$Group = factor(dat.plot$group,levels = c('AL','ICH'))
dat.plot$SampleID = dat.plot$orig.ident
dat.plot <- dat.plot[,.(N=.N),by=c("SampleID","Celltype","Group")]
sum.data <- dat.plot[,{.(NTotal = sum(.SD$N))},by=c("SampleID","Group")]
dat.plot = left_join(dat.plot,sum.data)
dat.plot$freq = dat.plot$N / dat.plot$NTotal *100
head(dat.plot)
compraisions.list=list( c('AL','ICH'))
g.colSet = list(group = c("AL" = mycolor[1],"ICH" =mycolor[2]))
p.list.perMcls <- lapply(unique(sort(dat.plot$Celltype)),function(mcls){
  
  text.size = 10
  text.angle = 45
  text.hjust = 1
  legend.position = "none"
  p <- ggboxplot(dat.plot[Celltype== mcls ,],x="Group",y="freq",
                 color = "Group", legend="none",title=mcls,
                 xlab="",ylab="frequency",
                 add = "jitter",outlier.shape=NA) +
    scale_color_manual(values=g.colSet$group) +
    scale_color_manual(values= c("AL" = mycolor[1],"ICH" =mycolor[2])) +
    scale_fill_manual(values= c("AL" = mycolor[1],"ICH" =mycolor[2])) +
    stat_compare_means(label="p.format",comparisons=compraisions.list) +  
    coord_cartesian(clip="off") +
    theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = text.size,color ="black"), 
          axis.text = element_text(size=text.size,color = "black"),
          axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), 
          legend.position = legend.position,
          legend.text = element_text(size= text.size),
          legend.title= element_text(size= text.size),
          strip.background = element_rect(color="black",size= 1, linetype="solid") 
    )
  return(p)
})
p.prop = wrap_plots(p.list.perMcls,ncol = 7)
p.prop
ggsave(p.prop,filename = "mouse.cellratio.pdf",width =23,height =6,units = "cm")


load('humansce.Rdata')
Idents(sce) <- 'celltype'
ann=c("Neutrophil", "monocyte"  , "T" ,"NK", "B"    )
sce$celltype <- factor(sce$celltype,levels = ann)
gene <- c("RGS2", "SOD2", "BASP1",
          "MS4A7", "CD86",  "VCAN",
          "CD3G","CD3D",  "CD3E", 
          "KLRD1", "GZMB","ADGRG1", 
          "CD79A", "MS4A1", "BANK1" )
jjDotPlot(object = sce,assay='RNA',slot='data',
          gene = gene,  id = 'celltype', 
          split.by.aesGroup =T, 
          xtree = F,ytree = F,
          rescale = T,rescale.min = -2,rescale.max =2,   midpoint = 0, dot.col = c('grey','white','#E18727FF'),
          point.geom = T,point.shape = 21,tile.geom = F)
graph2pdf(file='humandotplotann.pdf',width=7,height=6)
SCP::CellDimPlot(
  srt = sce, group.by = "celltype",
  reduction = "UMAP", bg_color = "white",
  pt.size = 2,pt.alpha = 1,raster = T,
  palcolor=mycolor,
  split.by = 'group',
  label = F,label.size =4,label_insitu = TRUE,label_repel = F,
  add_density = F, 
  theme_use = "theme_blank"
)
graph2pdf(file='humandim.pdf',width=8,height=6)

meta.tb = sce@meta.data
dat.plot = data.table(meta.tb)
dat.plot$Celltype = dat.plot$celltype
dat.plot$Group = factor(dat.plot$group,levels = c('HBP','MB'))
dat.plot$SampleID = dat.plot$orig.ident
dat.plot <- dat.plot[,.(N=.N),by=c("SampleID","Celltype","Group")]
sum.data <- dat.plot[,{.(NTotal = sum(.SD$N))},by=c("SampleID","Group")]
dat.plot = left_join(dat.plot,sum.data)
dat.plot$freq = dat.plot$N / dat.plot$NTotal *100
head(dat.plot)
compraisions.list=list( c('HBP','MB'))
g.colSet = list(group = c("HBP" = mycolor[1],"MB" =mycolor[2]))
p.list.perMcls <- lapply(unique(sort(dat.plot$Celltype)),function(mcls){
  
  text.size = 10
  text.angle = 45
  text.hjust = 1
  legend.position = "none"
  p <- ggboxplot(dat.plot[Celltype== mcls ,],x="Group",y="freq",
                 color = "Group", legend="none",title=mcls,
                 xlab="",ylab="frequency",
                 add = "jitter",outlier.shape=NA) +
    scale_color_manual(values=g.colSet$group) +
    scale_color_manual(values= c("HBP" = mycolor[1],"MB" =mycolor[2])) +
    scale_fill_manual(values= c("HBP" = mycolor[1],"MB" =mycolor[2])) +
    stat_compare_means(label="p.format",comparisons=compraisions.list) + 
    coord_cartesian(clip="off") +
    theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = text.size,color ="black"), 
          axis.text = element_text(size=text.size,color = "black"),
          axis.text.x = element_text(angle = text.angle, hjust = text.hjust ),
          legend.position = legend.position,
          legend.text = element_text(size= text.size),
          legend.title= element_text(size= text.size),
          strip.background = element_rect(color="black",size= 1, linetype="solid") 
    )
  return(p)
})
p.prop = wrap_plots(p.list.perMcls,ncol = 7)
p.prop
ggsave(p.prop,filename = "human.cellratio.pdf",width =18,height =6,units = "cm")

####step.3 FindMarkers#####
load('humansce.Rdata')
#load('mousesce.Rdata')
library(Seurat)
library(dplyr)
pbmc<-sce
Idents(pbmc) <- 'celltype'
pbmc$celltype.group <- paste(pbmc$celltype, pbmc$group, sep = "_")
pbmc$celltype <- Idents(pbmc)
Idents(pbmc) <- "celltype.group"
cellfordeg<-levels(pbmc$celltype)
data_forplot <- data.frame() 
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(pbmc, ident.1 = paste0(cellfordeg[i],"_MB"), ident.2 = paste0(cellfordeg[i],"_HBP"),assay = 'RNA',slot = 'counts',
                         logfc.threshold =0,min.pct = 0.2, verbose = T)
  CELLDEG$gene=rownames(CELLDEG)
  CELLDEG$cluster=cellfordeg[i]
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
  data_forplot=rbind(CELLDEG,data_forplot)
}
list.files()
ann=c("Neutrophil", "monocyte"  , "T" ,"NK", "B"    )
data_forplot$cluster <- factor(data_forplot$cluster,levels = ann)
library(scRNAtoolVis)
markerVolcano(markers =data_forplot,log2FC=0.1,
              topn =10,nforce=10,
              labelCol =mycolor
              
)
graph2pdf(file = 'volDEGs.pdf',height = 6,width =12)

pdf(file = 'DEGs.pdf',height = 5.5,width =12)
jjVolcano(diffData = data_forplot,
          log2FC.cutoff = 0.25,
          aesCol = c("#0072B5FF","#BC3C29FF"),
          topGeneN = 10,
          tile.col =c("#4DBBD5FF","#E64B35FF","#00A087FF","#3C5488FF", "#F39B7FFF"),
          size  =2,  fontface = 'italic',
          flip = F, 
          polar =F, 
)+ylim(-5,5)  
dev.off()
########step.4 get signature######
setwd('/Users/drning/Desktop/ICH/2.scRNA注释/human中性粒细胞')
library(nichenetr)
library(dplyr)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(AUCell)
library(ggplot2)
library(Seurat)
library(clusterProfiler)
library(IOBR)
load('humansce.Rdata')
Idents(sce) <- 'group'
CELLDEG <- FindMarkers(sce, ident.1 = 'MB', 
                       ident.2 = 'HBP',
                       logfc.threshold =0.25,min.pct = 0.25, verbose = T)
DEG=filter(CELLDEG,p_val_adj<0.05)
upgene=rownames(filter(DEG,avg_log2FC>0.25))
downgene=rownames(filter(DEG,avg_log2FC< -0.25))
signature=list(signature_up=upgene,
               signature_down=downgene)
save(signature,file = 'humansignature.Rdata')
#load('humansignature.Rdata')
#to mouse
signature_up <- na.omit(nichenetr::convert_human_to_mouse_symbols(symbols = signature[["signature_up"]]))
signature_down <- na.omit(nichenetr::convert_human_to_mouse_symbols(symbols = signature[["signature_down"]]))
signature=list(signature_up=signature_up,
               signature_down=signature_down)
##CCCF
sce <- irGSEA.score(object = sce, 
                             seeds = 123, ncores =30,
                             min.cells = 3, min.feature = 0,
                             custom = T, geneset = signature, 
                             species = "Homo sapiens", 
                             geneid = "symbol",
                             aucell.MaxRank = NULL, ucell.MaxRank = 3000,
                             kcdf = 'Gaussian')

load('neu.Rdata')
CellDimPlot(
  srt = sce, group.by = c("cell"),palcolor =c("#E8E8E8","#0072B5FF","#F47D2B"),
  pt.size = 2,
  bg_color = 'white',raster = T,
  pt.alpha = 1,
  reduction = "UMAP", theme_use = "theme_blank"
)
graph2pdf(file='humanumap1.pdf',height=4,width=4)

load('humanneufin.Rdata')
sce=humanneu
CellStatPlot(srt = sce, plot_type = c('rose'),force = T,
             palcolor =c("#0072B5FF","#F47D2B"),
             stat.by = "cell", group.by = "orig.ident",
             label = F)
graph2pdf(file='细胞比例扇图.pdf',height=5,width=5)
###鼠细胞比例
setwd('/Users/drning/Desktop/ICH/2.scRNA注释/mouse中性粒细胞')
load('neu.Rdata')
CellDimPlot(
  srt = sce, group.by = c("cell"),palcolor =c("#E8E8E8","#F47D2B","#0072B5FF"),
  pt.size = 2,
  bg_color = 'white',raster = T,
  pt.alpha = 1,
  reduction = "UMAP", theme_use = "theme_blank"
)
graph2pdf(file='mouseumap1.pdf',height=4,width=4)

load('mouseneufin.Rdata')
sce=mouseneu
CellStatPlot(srt = sce, plot_type = c('rose'),force = T,
             palcolor =c("#0072B5FF","#F47D2B"),
             stat.by = "cell", group.by = "orig.ident",
             label = F)
graph2pdf(file='细胞比例扇图.pdf',height=5,width=5)
###step.5 Pseudo-time######
#seurat2h5ad
seurat2h5ad<- function(data=sc){
  library(Seurat)
  library(SeuratDisk)
  seurat_object=data
  seu = DietSeurat(
    seurat_object,
    counts = TRUE, 
    data = TRUE, 
    scale.data = FALSE, 
    features = rownames(seurat_object), 
    assays = "RNA",
    dimreducs = c("pca","umap"),
    graphs = c("RNA_nn", "RNA_snn"), 
    misc = TRUE
  )
  i <- sapply(seu@meta.data, is.factor)
  seu@meta.data[i] <- lapply(seu@meta.data[i], as.character)
  SaveH5Seurat(seu, filename = "./SeuratDisk.h5seurat", overwrite = TRUE)
  Convert("./SeuratDisk.h5seurat", "./SeuratDisk.h5ad", assay="RNA", overwrite = TRUE)
}
sc=sce
scobj <- CreateSeuratObject(counts = sc@assays[["RNA"]]@counts,meta.data = sc@meta.data)
seurat2h5ad(data=scobj)
#run in python
import importlib_metadata
import sys
sys.modules['importlib.metadata'] = importlib_metadata
import sctour as sct
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import skmisc
import os
import anndata as ad
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
sc.settings.set_figure_params(dpi=50, facecolor="white")
new_path = "/Users/drning/Desktop/ICH/2.scRNA注释/"
os.chdir(new_path)
adata = sc.read('SeuratDisk.h5ad') 
adata
adata.X = adata.X.astype(int)
adata.raw = adata.copy()
sc.pp.calculate_qc_metrics(adata,  percent_top=None, log1p=False, inplace=True)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=1000,subset=True)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata,color="celltype",size=2)
adata = adata.raw.to_adata()
adata.X = adata.X.astype('float32')
tnode = sct.train.Trainer(adata, loss_mode='nb', alpha_recon_lec=0.5, alpha_recon_lode=0.5)
tnode.train()
adata.obs['ptime'] = tnode.get_time()
mix_zs, zs, pred_zs = tnode.get_latentsp(alpha_z=0.5, alpha_predz=0.5)
adata.obsm['X_TNODE'] = mix_zs
adata.obsm['X_VF'] = tnode.get_vector_field(adata.obs['ptime'].values, adata.obsm['X_TNODE'])
adata = adata[np.argsort(adata.obs['ptime'].values), :]
sc.pp.neighbors(adata, use_rep='X_TNODE', n_neighbors=15)
sc.tl.umap(adata)
sc.pl.umap(adata,color="celltype",size=2)
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(11, 11))
sc.pl.umap(adata, color='celltype', ax=axs[0, 0], legend_loc='on data', show=False, frameon=False)
sc.pl.umap(adata, color='ptime', ax=axs[1, 0], show=False, frameon=False)
sct.vf.plot_vector_field(adata, zs_key='X_TNODE', vf_key='X_VF', use_rep_neigh='X_TNODE', 
                         color='ptime', show=False, ax=axs[1, 1], legend_loc='none', frameon=False, size=100, alpha=0.2,
                         save='sctour.pdf')
plt.show()
sc_embedding = pd.DataFrame(adata.obsm["X_umap"],columns=['umap_1', 'umap_2'],index=adata.obs.index )
sc_embedding.to_csv('sc_embedding.csv',sep='\t')
ptime_df = adata.obs[['ptime']]
ptime_df.to_csv('ptime.csv', index=True)
#run in R
sc_embedding <- read.table('sc_embedding.csv',header = T,row.names = 1,sep='\t')
head(sc_embedding)
sc[["umap"]] <- CreateDimReducObject(as.matrix(sc_embedding),key = "umap_")
ptime=read.csv('ptime.csv')
rownames(ptime)=ptime[,1]
ptime=ptime[colnames(sc),]
sc$ptime=ptime$ptime
sc$celltype=as.character(sc$celltype)

sc <- RunDynamicFeatures(srt = sc, lineages = c("ptime"), 
                          n_candidates = 2000,slot = "data")
load('gene.Rdata')
ht <- DynamicHeatmap(
  srt = sc, lineages = c("ptime"),,slot = "data",assay = 'magic',
  use_fitted = TRUE, n_split = 3, reverse_ht = "ptime",
  species = "Homo_sapiens", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "magma", cell_annotation = "cell",exp_method='zscore',topTerm=10,
  features_label=gene,
  separate_annotation_palette = c("Paired", "Set1"),
  pseudotime_label = 25,
  height = 5, width = 2
)
print(ht$plot)
graph2pdf(file="heatmap.pdf",width=15,height=8)
save(ht,file = 'ht.Rdata')
###step.6 SCENIC######
#prepare input data
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
easy_scenic <- function(scdata=sce,
                        celltype=sce$siggroup,
                        org="hgnc",   
                        nCores=35,
                        dbDir="/home/ningjingyuan/SCENICref"
){
  library(plyr)
  library(permute)
  library(data.table)
  library(SCopeLoomR)
  library(SCENIC)
  sce <- scdata
  sce$CellState <- celltype
  cell.info <- sce@meta.data
  head(cell.info)
  dge <- as.matrix(sce@assays$RNA@counts)
  if (org == "hgnc") {
    mydbs <- c("hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather","hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
  }
  if (org == "mgi") {
    mydbs <- c("mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather", "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
  }
  names(mydbs) <- c("10kb","500bp")
  scenicOptions <- initializeScenic(org=org,  
                                    nCores=nCores,
                                    dbDir=dbDir, 
                                    dbs = mydbs,
                                    datasetTitle = "os")
  genesKept <- geneFiltering(
    exprMat=dge, 
    scenicOptions=scenicOptions,
    minCountsPerGene = 1,
    minSamples = 20
  )
  length(genesKept)
  ## [1] 18250
  dge <- dge[genesKept, ]
  dim(dge)
  source("/home/ningjingyuan/codesource/scenic/AvgN.R")
  avg20.rep1 <- AvgN(dge, cell.info, seed=1)
  #avg20.rep2 <- AvgN(dge, cell.info, seed=2)
  # avg20.rep3 <- AvgN(dge, cell.info, seed=3)
  if(!dir.exists("output")) {
    dir.create("output")
  }
  saveRDS(cell.info, "output/s1_cell.info.rds")
  source("/home/ningjingyuan/codesource/scenic/add_cellAnnotation.R")
  saveLoom <- function(exprMat, output){
    ## prepare cell cluster info
    cellInfo <- data.frame(
      row.names = colnames(exprMat),
      #cellType = sapply(strsplit(colnames(exprMat), split = "\\."), function(x) x[1])
      cellType = sub("\\.[0-9]+", "", colnames(exprMat))
    )
    ## 存储loom文件
    loom <- build_loom(output, dgem=exprMat)
    loom <- add_cellAnnotation(loom, cellInfo)
    close_loom(loom)
  }
  
  saveLoom(avg20.rep1, "output/s1_avg20_rep1.loom")
  #saveLoom(avg20.rep2, "output/s1_avg20_rep2.loom")
  # saveLoom(avg20.rep3, "output/s1_avg20_rep3.loom")
  loom <- build_loom("output/s1_exprMat.loom", dgem=dge)
  loom <- add_cellAnnotation(loom, cell.info)
  close_loom(loom)
  
}
easy_scenic(scdata=sc,
            celltype=sc$celltype,
            org="hgnc",   
            nCores=35,
            dbDir="/home/ningjingyuan/SCENICref")
#run in linux
conda activate scenic_protocol 
cd /scenic
output_path="/scenic/output"
nohup bash "$output_path/s2_runPySCENIChuman.sh" avg20_rep1 > run.log 2>&1 &
f_loom_path_scenic="$output_path/s2_avg20_rep1.pyscenic.loom"
ctx_output="$output_path/s2_avg20_rep1.reg.tsv"
sample_name="$output_path/s3_avg20_rep1"
threads=30
min_regulon_size=10
nohup python /home/ningjingyuan/codesource/scenic/s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size &


#plot
load('humanneufin.Rdata')
sce=humanneu
CellDimPlot(
  srt = sce, group.by = c("cell"),palcolor =c("#E8E8E8","#F47D2B","#0072B5FF"),
  pt.size = 2,
  bg_color = 'white',raster = T,
  pt.alpha = 1,
  reduction = "UMAP", theme_use = "theme_blank"
)
data=as.data.frame(sce@assays[["scenic"]]@data)
data['ptime',]=sce$ptime
sceniccor=easy_corbatch(data=data,
                        gene='ptime', 
                        needplot=F, 
                        corFilter=0.3,
                        pFilter=0.05
)
###step.7 ML##########
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(circlize)
Cindex_mat <- read.table("/Users/drning/Desktop/ICH/3ML/Results/AUC_mat.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
avg_Cindex <- sort(rowMeans(Cindex_mat[, 2:3]), decreasing = TRUE)
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]  
avg_Cindex <- round(avg_Cindex, 3) 
CohortCol <- c("#90D5E4", "#89C75F", "#F37B7D")
names(CohortCol) <- colnames(Cindex_mat)
color_fun <- colorRamp2(c(0.7, 0.85, 1), c("#3BBCA8", "#FFFFFF", "#D24B27"))
color_fun <- colorRamp2(c(0.7, 0.85, 1), c("#3BBCA8", "#FFFFFF", "#D24B27"))
col_ha <- columnAnnotation(Cohort = colnames(Cindex_mat), col = list(Cohort = CohortCol), show_annotation_name = FALSE)
row_ha <- rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                           gp = gpar(fill = "steelblue", col = NA),
                                           add_numbers = TRUE, numbers_offset = unit(-10, "mm"),
                                           axis_param = list(labels_rot = 0),
                                           numbers_gp = gpar(fontsize = 9, col = "white"),
                                           width = unit(3, "cm")),
                        show_annotation_name = FALSE)
# 绘制热图
heatmap <- Heatmap(as.matrix(Cindex_mat), name = "AUC",
                   right_annotation = row_ha, 
                   top_annotation = col_ha,
                   col = color_fun, 
                   rect_gp = gpar(col = "black", lwd = 1),  
                   cluster_columns = FALSE, 
                   cluster_rows = FALSE,  
                   show_column_names = FALSE, 
                   show_row_names = TRUE,
                   row_names_side = "left",
                   width = unit(1 * ncol(Cindex_mat) + 2, "cm"),
                   height = unit(0.5 * nrow(Cindex_mat), "cm"),
                   column_split = c(rep(1,1),rep(2,2)), 
                   column_title = NULL,
                   cell_fun = function(j, i, x, y, w, h, col) {  
                     grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                               x, y, gp = gpar(fontsize = 10))
                   }
)
pdf("/Users/drning/Desktop/ICH/3ML/Figures/heatmap.pdf", width = 8, height = 30)
draw(heatmap)
dev.off()

model=readRDS('/Users/drning/Desktop/ICH/3ML/Results/model.rds')
res=as.data.frame(coef(model[["Lasso+Stepglm[forward]"]]))
colnames(res)='coef'
res$gene=rownames(res)
res=res[-1,]
res=arrange(res,desc(coef))
data=res
data$gene = factor(data$gene,
                     levels = unique(data$gene))
Color = easy_color(palette_num = 45,n = 7)
ggplot(data, aes(x = gene, y = coef,colour  = gene )) +
  geom_point(aes(color = gene), size = 1.5) +
  geom_point(aes(color = gene), size = 3, shape = 21, fill = NA) +
  geom_segment(aes(xend = gene, y = 0, yend = coef, color = gene), linewidth = 1) +
  scale_color_manual(values = Color) +
  theme_classic() +
  theme(legend.position = "none",  
        axis.text.x = element_text(size = 9,
                                   colour = "black",
                                   vjust = 0.5,
                                   hjust = 0.5,
                                   angle = 90),
        axis.line = element_line(color = "black"),
  ) +
  labs(x = "", y = "Coef")
graph2pdf(file='coef.pdf',width=1.5,height=1.5)

#### Step 8: Experimental Plotting ####
###发病曲线###
sur=read.csv('/发病曲线.csv')
plotdata=data.frame(row.names = rownames(sur),
                    fustat=sur$status,
                    futime=sur$time,
                    group=sur$group)
head(plotdata)
plotdata=filter(plotdata,group %in% c('Dihydroergotamine','Vehicle'))
surdata <-plotdata
diff <- survdiff(Surv(futime, fustat) ~ group, data = surdata)
length <- length(levels(factor(surdata$group)))
pValue <- 1 - pchisq(diff$chisq, df = length - 1)
print(paste0('p=',pValue))
fit <- survfit(Surv(futime, fustat) ~ group, data = surdata)

#Kaplan-Meier生存曲线
ggsurvplot(
  fit,
  data = surdata,
  pval = pValue,
  pval.size = 3.5,
  legend.labs = c("Dihydroergotamine","Vehicle"),
  palette = c( "#1F78B4","#A6CEE3"),
  break.time.by = 5,
  conf.int = F,
  risk.table =F,
  risk.table.title = 'group',
  ggtheme = theme_classic() + theme(panel.grid = element_blank()), 
  main = "Survival curve"
)
graph2pdf(file='发病曲线.pdf',width=3,height=3)
###脑出血量
library(ggplot2)
library(tidyverse)
library(plyr)
library(ggbeeswarm)
library(ggsignif)
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
data=read.csv('脑出血量.csv')
head(data)
data=data[,c(1,2)]
dt_long <- gather(data,group,value)
head(dt_long)
str(dt_long)
dt_long$group <- factor(dt_long$group,levels = c('vehicle','Dihydroergotamine'))
dt_sum <- data_summary(dt_long, varname="value",groupnames="group")
head(dt_sum)
# draw
ggplot() + 
  geom_bar(dt_sum,mapping = aes(x=group, y=value, fill=group),
           stat="identity",width=.6,alpha = 0.8,position=position_dodge()) +
  geom_errorbar(dt_sum,mapping = aes(x=group, y=value,ymin=value-sd, ymax=value+sd), 
                width=.3,position=position_dodge(.8)) +
  geom_beeswarm(dt_long,mapping = aes(x=group, y=value, fill=group),
                shape = 21,color = 'black',size = 2,cex = 3,stroke = 0.6)+
  geom_signif(dt_long,mapping = aes(x=group, y=value),
              comparisons = list(c("vehicle","Dihydroergotamine")),
              test = "t.test",
              step_increase = 0.1,
              map_signif_level = F)+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4"))+
  labs(x = "group", y = "Cerebral hemorrhage volume(ul)")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14, angle = 25, hjust = 1),
        axis.title.y=element_text(color="black",size=14))
ggsave('出血体积.pdf',width = 4,height =4) 

data=read.csv('炎症细胞浸润.csv')
head(data)
data=data[,c(1,2)]
dt_long <- gather(data,group,value)
head(dt_long)
str(dt_long)
dt_long$group <- factor(dt_long$group,levels = c('Vehicle','Dihydroergotamine'))
dt_sum <- data_summary(dt_long, varname="value",groupnames="group")
head(dt_sum)
# draw
ggplot() + 
  geom_bar(dt_sum,mapping = aes(x=group, y=value, fill=group),
           stat="identity",width=.6,alpha = 0.8,position=position_dodge()) +
  geom_errorbar(dt_sum,mapping = aes(x=group, y=value,ymin=value-sd, ymax=value+sd), 
                width=.3,position=position_dodge(.8)) +
  geom_beeswarm(dt_long,mapping = aes(x=group, y=value, fill=group),
                shape = 21,color = 'black',size = 2,cex = 3,stroke = 0.6)+
  geom_signif(dt_long,mapping = aes(x=group, y=value),
              comparisons = list(c("Vehicle","Dihydroergotamine")),
              test = "t.test",
              step_increase = 0.1,
              map_signif_level = F)+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4"))+
  labs(x = "group", y = "Inflammatory cell infiltration")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14, angle = 25, hjust = 1),
        axis.title.y=element_text(color="black",size=14))
ggsave('炎症细胞浸润量.pdf',width = 4,height =4) 


data=read.csv('氧化应激.csv')
head(data)
data=data[,c(1,2)]
dt_long <- gather(data,group,value)
head(dt_long)
str(dt_long)
dt_long$group <- factor(dt_long$group,levels = c('Vehicle','Dihydroergotamine'))
dt_sum <- data_summary(dt_long, varname="value",groupnames="group")
head(dt_sum)
# draw
ggplot() + 
  geom_bar(dt_sum,mapping = aes(x=group, y=value, fill=group),
           stat="identity",width=.6,alpha = 0.8,position=position_dodge()) +
  geom_errorbar(dt_sum,mapping = aes(x=group, y=value,ymin=value-sd, ymax=value+sd), 
                width=.3,position=position_dodge(.8)) +
  geom_beeswarm(dt_long,mapping = aes(x=group, y=value, fill=group),
                shape = 21,color = 'black',size = 2,cex = 3,stroke = 0.6)+
  geom_signif(dt_long,mapping = aes(x=group, y=value),
              comparisons = list(c("Vehicle","Dihydroergotamine")),
              test = "t.test",
              step_increase = 0.1,
               map_signif_level = F)+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4"))+
  labs(x = "group", y = "oxidative stress")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14, angle = 25, hjust = 1),
        axis.title.y=element_text(color="black",size=14))
ggsave('oxidative stress.pdf',width = 4,height =4) 


data=read.csv('tune.csv')
head(data)
data=data[,c(1,2)]
dt_long <- gather(data,group,value)
head(dt_long)
str(dt_long)
dt_long$group <- factor(dt_long$group,levels = c('Vehicle','Dihydroergotamine'))
dt_sum <- data_summary(dt_long, varname="value",groupnames="group")
head(dt_sum)
# draw
ggplot() + 
  geom_bar(dt_sum,mapping = aes(x=group, y=value, fill=group),
           stat="identity",width=.6,alpha = 0.8,position=position_dodge()) +
  geom_errorbar(dt_sum,mapping = aes(x=group, y=value,ymin=value-sd, ymax=value+sd), 
                width=.3,position=position_dodge(.8)) +
  geom_beeswarm(dt_long,mapping = aes(x=group, y=value, fill=group),
                shape = 21,color = 'black',size = 2,cex = 3,stroke = 0.6)+
  geom_signif(dt_long,mapping = aes(x=group, y=value),
              comparisons = list(c("Vehicle","Dihydroergotamine")),
              test = "t.test",
              step_increase = 0.1,
               map_signif_level = F)+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4"))+
  labs(x = "group", y = "tunel")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14, angle = 25, hjust = 1),
        axis.title.y=element_text(color="black",size=14))
ggsave('tunel.pdf',width = 4,height =4) 


data=read.csv('/中性粒细胞.csv')
head(data)
data=data[,c(1,2)]
dt_long <- gather(data,group,value)
head(dt_long)
str(dt_long)
dt_long$group <- factor(dt_long$group,levels = c('Vehicle','Dihydroergotamine'))
dt_sum <- data_summary(dt_long, varname="value",groupnames="group")
head(dt_sum)
# draw
ggplot() + 
  geom_bar(dt_sum,mapping = aes(x=group, y=value, fill=group),
           stat="identity",width=.6,alpha = 0.8,position=position_dodge()) +
  geom_errorbar(dt_sum,mapping = aes(x=group, y=value,ymin=value-sd, ymax=value+sd), 
                width=.3,position=position_dodge(.8)) +
  geom_beeswarm(dt_long,mapping = aes(x=group, y=value, fill=group),
                shape = 21,color = 'black',size = 2,cex = 3,stroke = 0.6)+
  geom_signif(dt_long,mapping = aes(x=group, y=value),
              comparisons = list(c("Vehicle","Dihydroergotamine")),
              test = "t.test",
              step_increase = 0.1,
               map_signif_level = F)+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4"))+
  labs(x = "group", y = "neutrophil")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14, angle = 25, hjust = 1),
        axis.title.y=element_text(color="black",size=14))
ggsave('neutrophil.pdf',width = 4,height =4) 
###出血灶###
library(ggpubr)
library(dplyr)
library(rstatix)
data <- read.csv('出血灶数量.csv')
stat.test <- data %>%
  group_by(avg) %>%
  t_test(number ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "group")
custom_colors <- c("Vehicle" = "#A6CEE3", "Dihydroergotamine" =  "#1F78B4")
ggbarplot(
  data, x = "group", y = "number", fill = 'group',
  add = c("mean_sd", "jitter"), facet = c("avg")
) +
  stat_pvalue_manual(stat.test) +
  scale_fill_manual(values = custom_colors) + 
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_line(size = 0.6),
    axis.ticks.y = element_line(size = 0.6),
    axis.title = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(vjust = -0.3, size = 12),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_blank()
  )
ggsave('出血灶.pdf',width = 4,height =7 ) 
####sbp####
library(dplyr)
library(ggplot2)
library(plotrix)
library(tidyverse)
library(reshape2)
data<- read.csv('血压曲线.csv')
data=filter(data,treat %in% c('Dihydroergotamine','Vehicle'))
data1<-melt(data,id=c('treat','Day'))
topbar<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottombar<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
data1$treat=factor(data1$treat,levels = c(c('Vehicle','Dihydroergotamine')))
ggplot(data1,aes(Day,value,color=treat))+
  stat_summary(geom='line',fun='mean',cex=1)+
  stat_summary(geom='errorbar',
               fun.min=bottombar,fun.max=topbar,
               width=1,cex=0.8,aes(color=treat))+
  stat_summary(geom='point',fun='mean',aes(fill=treat),
               size=2,pch=21,color='black')+
  theme_classic(base_size=15)+
  scale_fill_manual(values=c("#A6CEE3","#1F78B4"))+
  scale_color_manual(values=c("#A6CEE3","#1F78B4"))+
  labs(x = "Time (days)", y = 'SBP (mmHg)')
ggsave('SBP.pdf',width = 4.7,height =3 ) 

data=read.csv('/AL.ICH.Neu.csv')
bar.plot1(data = data,group = 'Group',feature = 'Percentage',
          group.levels = c('AL','ICH'),compar.groups =list(c('AL','ICH') ),point.size = 2.5,
          x.label = 'Group',y.label = 'NXPE3+ neutrophil(%)' )
ggsave('AL.ICH.Neu.pdf',width = 2.5,height =3 ) 
data=read.csv('Veh.Dih.Neu.csv')
bar.plot1(data = data,group = 'Group',feature = 'Percentage',
          group.levels = c('Veh','Dih'),compar.groups =list( c('Veh','Dih') ),point.size = 2.5,
          x.label = 'Group',y.label = 'NXPE3+ neutrophil(%)' )
ggsave('Veh.Dih.Neu.pdf',width = 2.5,height =3 ) 


data=read.csv('Al.ICH.WB.csv')
bar.plot1(data = data,group = 'Group',feature = 'plot.numer',
          group.levels = c('AL','ICH'),compar.groups =list(c('AL','ICH') ),point.size = 2.5,
          x.label = 'Group',y.label = 'NXPE3(fold of AL)' )
ggsave('AL.ICH.WB.pdf',width = 2.5,height =3 ) 
data=read.csv('Veh.Dih.WB.csv')
bar.plot1(data = data,group = 'Group',feature = 'plot.numer',
          group.levels = c('Veh','Dih'),compar.groups =list( c('Veh','Dih') ),point.size = 2.5,
          x.label = 'Group',y.label = 'NXPE3(fold of Veh)' )
ggsave('Veh.Dih.WB.pdf',width = 2.5,height =3 ) 



















































