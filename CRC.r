library(ggplot2)
library(limma) 
library(glmnet)
library(survival)
library(dplyr)
library(ggpubr)
library(survminer)
library(rtracklayer)
library(clusterProfiler)
library(enrichplot)
library(Hmisc)
library(pheatmap)
library(ggVolcano)
library(RColorBrewer)
library(eoffice)
library(colorspace) 

mRNA_TPM<-GetRowNames(read.csv("mRNA_TPM.csv",header = T))
colnames(mRNA_TPM)<-gsub("X","",colnames(mRNA_TPM))
LNC_TPM<-GetRowNames(read.csv("LNC_TPM.csv",header = T))
colnames(LNC_TPM)<-gsub("X","",colnames(LNC_TPM))
TCGA_TPM<-GetRowNames(read.csv("TCGA_TPM.csv",header = T))

library(clusterProfiler)
library(enrichplot)
tmp=read.table("GeneList.txt",sep="\t",header = T,stringsAsFactors = F)
geneSet<-unique(tmp$Category)
  head(tmp)
 write.gmt <- function(geneSet=geneSet,gmt_file='immPort.gmt'){
 sink( gmt_file )
 for (i in 1:length(geneSet)){
 cat(geneSet[i])
 cat('\tNA\t')
 cat(paste(tmp$Symbol[which(tmp$term==geneSet[i])],collapse = '\t'))
 cat('\n') }
 sink()
 }
 write.gmt(geneSet,'immPort.gmt')
 ###新版生产了 immune2
S<-read.gmt(gmtfile = "all_gene_lists.gmt") 
 ###旧版生产了 immune  后期复现产生了immune3
S<-tmp[,c("Category","Symbol")]
head(S)
colnames(S)<-c("term","gene")
length(unique(S$gene))
table(S$term)
res<-cor(t(LNC_TPM),t(mRNA_TPM),method="spearman")
immune<-IMCUTOFF(res)
write.csv(immune,"immune3.csv")
immune_Gene<- tmp$Symbol
length(immune_Gene)





immune<-read.csv("./immune3.csv",header = T,row.names = 1)
Get_label <- function(EXP) {
    Temp <- as.matrix(rep(NA, times = (dim(EXP)[2])))
    ##### 自定义####
    Temp[grep("N", rownames(t(EXP)))] <- 1
    Temp[grep("T", rownames(t(EXP)))] <- 2
    colnames(Temp) <- "label"
    Temp2 <- rownames(t(EXP))
    Clinic <- as.data.frame(cbind(Temp2, Temp))
    return(Clinic)
}
Get_label2 <- function(EXP) {
    Temp <- as.matrix(rep(NA, times = (dim(EXP)[2])))
    ##### 自定义####
    substr(colnames(EXP),14,16)!='11A'
    Temp[grep("11A", rownames(t(EXP)))] <- 1
    Temp[grep(".01", rownames(t(EXP)))] <- 2
    colnames(Temp) <- "label"
    Temp2 <- rownames(t(EXP))
    Clinic <- as.data.frame(cbind(Temp2, Temp))
   # Clinic$Temp<-as.numeric(Clinic$Temp)
    return(Clinic)
}
foldChange<-function(KKK,i,K){
   fpkm<-KKK
#tpms <- apply(fpkm,2,fpkmToTpm)
        if(sum(is.na(K))!=0){
        K<-na.omit(K)
        fpkm<-fpkm[,K[,1]] 
        K<-K[,2]
                        }
group_list<-K

exprSet <- fpkm
exprSet=normalizeBetweenArrays(exprSet)
if(i==0){#print(head(exprSet))
        }#查看数据是否要转换
if(i==1){exprSet <- log2(exprSet+1)}
#差异分析：
dat <- exprSet
design=model.matrix(~factor( group_list ))
colnames(design) <- c("Control", "Treatment")
fit=lmFit(dat,design)
contrast_matrix <- makeContrasts(Treatment - Control, levels=design)
fit2 <- contrasts.fit(fit,  contrast_matrix)
fit=eBayes(fit)
#options(digits = 4)
#topTable(fit,coef=2,adjust='BH')
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
    deg$row<-rownames(deg)
    return(deg)
}
L2<-foldChange(LNC_TPM,1,Get_label(LNC_TPM)[[2]])
#C2<-foldChange(circRNA,1,Get_label(LncRNA)[[2]])
head(L2)
L3<-foldChange(mRNA_TPM,1,Get_label(LNC_TPM)[[2]])
#C2<-foldChange(circRNA,1,Get_label(LncRNA)[[2]])
head(L3)
L4<-foldChange(TCGA_TPM,0,Get_label2(TCGA_TPM))
#C2<-foldChange(circRNA,1,Get_label(LncRNA)[[2]])
head(L4)

filterIMDE <- function(immune, L2,L, TCGA_LNC2,IMCUTOFF, DECUTOFF,L4) {
    RESULT<-list()
    temp1 <- as.character(unique(immune[, 1]))
    immuneGENE <- NULL
    for (i in 1:length(temp1)) {
        temp <- immune[which(immune[, 1] == temp1[i]), ]
        temp2 <- as.character(unique(temp[which(abs(temp$enrichmentScore) > IMCUTOFF[1] & 
            temp$p.adjust < IMCUTOFF[2]), 1]))
        immuneGENE <- c(immuneGENE, temp2)
    } 
    {print("immuneGENE Number:")
    print(length(immuneGENE))
    LLL<-rownames(L2)[which(abs(L2$logFC)>=DECUTOFF[1]&L2$adj.P.Val<DECUTOFF[2])]
    deGENE<-intersect(LLL,rownames(L))
    print("deGENE Number:")
    print(length(deGENE))
    deimGENE<-intersect(deGENE,immuneGENE)
    print("deimGENE Number:")
    print(length(deimGENE))
    FdeimGENE<-intersect(rownames(TCGA_LNC2),deimGENE)
    print("TCGAdeimGENE Number:")
    print(length(FdeimGENE))}
    LLL4<-rownames(L4)[which(abs(L4$logFC)>=DECUTOFF[1]&L4$adj.P.Val<DECUTOFF[2])]
    LLI<-intersect(LLL4,FdeimGENE)
    print(length(LLI))
    RESULT[[1]]<-FdeimGENE
    RESULT[[2]]<-L[match(FdeimGENE,rownames(L)),]
    RESULT[[3]]<-TCGA_LNC2[match(FdeimGENE,rownames(TCGA_LNC2)),]
    RESULT[[4]]<-L[match(deimGENE,rownames(L)),]
    RESULT[[5]]<-immuneGENE
    RESULT[[6]]<-deimGENE
    
    names(RESULT)<-c("FDIG","FDIG_LNC_TMP","TCGA_FDIG_LNC_TMP","DIG_LNC_TMP")
    return(RESULT)
    }


LNC<-filterIMDE(immune,L2,LNC_TPM,TCGA_TPM,c(0.5,0.01),c(1.0,0.01),L4)
names(LNC)
indata <- LNC[[3]][, substr(colnames(LNC$TCGA_FDIG_LNC_TMP), 14, 16) != "11A"]%>%scale()
SV$match_id <- gsub("-", ".", SV[, 1])



temp<-immune[order(immune[,2],decreasing = T),]
deimGENE<-rownames(LNC[[4]])
FF<-temp[as.character(temp[,1])%in%deimGENE,]
FJ<-cbind(as.character(FF[,1]),gsub("[[:digit:]]","",rownames(FF)))
FJ1<-as.matrix(table(as.data.frame(FJ)))
colnames(FJ1)[3]<-'BCR Signaling Pathway'
colnames(FJ1)[9]<-'Natural Killer_Cell_Cytotoxicity'
colnames(FJ1)[10]<-'TCR Signaling Pathway'
colnames(FJ1)[10]<-'TCR Signaling Pathway'
library("circlize")
#pdf("./plot/circos_plot3.pdf",width =55 ,height = 45)
#png("./plot/circos_plot3.png",width =6000 ,height = 4000)
i=45
set.seed(i)
#win.graph(width=20, height=20,pointsize=12)
 # plot chord siagram
cir_P<-as.ggplot(function(){
    options(repr.plot.width=8, repr.plot.height=8)
      # par(cex = 4, mar = c(0, 0, 0, 0))
             chordDiagram(FJ1,
             directional = 1,
              annotationTrack = c("grid"),
              annotationTrackHeight = c(0.1, 0.1),
              diffHeight  =-0.04,
              preAllocateTracks = 1,
             order = c(rownames(FJ1),colnames(FJ1)))
circos.track(track.index = 1, panel.fun = function(x, y) {
   circos.text(CELL_META$xcenter, CELL_META$ylim[1],CELL_META$sector.index,
               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
 }, bg.border = NA)
})

cir_P



TDATA<-as.data.frame(cbind(t(LNC[[4]]),Get_label(LNC[[4]])[2]))
TDATA[dim(TDATA)[2]]<-as.integer(unlist(TDATA[dim(TDATA)[2]]))
TDATA[dim(TDATA)[2]][which(TDATA[dim(TDATA)[2]]==2),]<-0
TDATA[dim(TDATA)[2]]<-as.integer(unlist(TDATA[dim(TDATA)[2]]))
head(TDATA)
NMS<-colnames(TDATA)[1:dim(TDATA)[2]-1]
rownames(TDATA)<-NULL
colnames(TDATA)<-paste("P",1:dim(TDATA)[2],sep = "")
K<-dim(TDATA)[2]
library(EFS)
efs <- ensemble_fs(data = log2(TDATA+1),classnumber = K,runs=10,
                   selection = c(TRUE,TRUE, TRUE,TRUE, TRUE, TRUE, FALSE, FALSE))
				   
EFSM2<-as.matrix(apply(efs,2,sum)[order(apply(efs,2,sum),decreasing = FALSE)])
colnames(efs)<-NMS

EFSM<-as.matrix(apply(efs,2,sum)[order(apply(efs,2,sum),decreasing = FALSE)])

EFSM2<-as.matrix(EFSM2[which(EFSM>0.46),])
EFSM<-as.matrix(EFSM[which(EFSM>0.46),])

KK<-matrix(0,nrow=length(rownames(EFSM)),ncol=length(NE[,1]))
rownames(KK)<-rev(rownames(EFSM))
colnames(KK)<-NE[,1]
PP<-KK
head(PP)
for(i in 1:length(rownames(EFSM)))
    {temp<-Q[which(Q[,1]==rev(rownames(EFSM))[i]),]
    KK[i,match(gsub("[[:digit:]]","",rownames(temp)),colnames(KK))]<-as.numeric(temp[,2])
    PP[i,match(gsub("[[:digit:]]","",rownames(temp)),colnames(KK))]<-as.numeric(temp[,4])
}
library("ggplotify")
library("corrplot")

P_cor<-as.ggplot(function()print(corrplot(t(KK[,1:15]), p.mat = t(PP[,1:15]),  method ="circle",sig.level = 0.05, insig = "p-value",
        col=col2(100),tl.col="black")))+coord_flip()
		
		
library(factoextra)
library(survminer)
library(ggplot2)
library(survivalROC)
library(survival)
library(ggplotify)

indata<-scale(LNC[[4]][rev(row.names(EFSM)),grep("T",colnames(LNC[[4]]))])
as.ggplot(fviz_nbclust(t(indata), kmeans, method = "wss") + geom_vline(xintercept = 2, linetype = 2))
km_result <- kmeans(t(indata), 2, nstart = 100)
table(km_result$cluster)
res<-km_result$cluster


source("/mnt//RA//whn//PAGE.r")
clic=data.frame("Pid"=names(res),"label"=c(res))
REO<-REOF1(mRNA_TPM[co_mRNA,],clic,7000,1.0e-5,15)



library(colorspace)
options(repr.plot.width=8, repr.plot.height=8)
i=15
theme_set(theme_gray(base_size = i)) 
theme1<-theme(
    plot.title = element_text(size = i),  # 主图标题
    plot.subtitle = element_text(size = i),
    plot.caption = element_text(size = i),
    axis.text = element_text(size = i),    # 主图坐标轴标签
    axis.title = element_text(size = i),
    legend.title = element_text(size = i),
    legend.text = element_text(size = i),
  )


indata<-scale(LNC[[4]][FF,grep("T",colnames(LNC[[4]]))])
F3.1<-as.ggplot(fviz_nbclust(t(indata), kmeans, method = "wss") + 
                geom_vline(xintercept = 2, linetype = 2))+theme(
    plot.title = element_text(size = i),  # 主图标题
    axis.text.x = element_text(size = 0),    # 主图坐标轴标签
    axis.text.y = element_text(size = 0),    # 主图坐标轴标签
    )
km_result <- kmeans(t(indata), 2, nstart = 100)
table(km_result$cluster)
res<-km_result$cluster
F3.1
#hcl_palettes(plot = TRUE)

index<-order(res)
ann_colors <- list(Cluster = c("1" = cols[5],"2" = cols[6]))
annotation_col = data.frame(
 Cluster = factor(res[index])
#TumorStage =factor(na.omit(FSV)$tumor_stage2[index])
)
#annotation_col
rownames(annotation_col) = colnames(indata[,index])
dist.r<-dist(t(indata[,index]),method="euclidean")
F3.2<-as.ggplot(function(){pheatmap(as.matrix(dist.r),
        cluster_row = FALSE,
         cluster_cols = FALSE,
                                     fontsize=10,
          show_colnames=F,
          show_rownames=F,
        legend = FALSE,
        annotation_legend = FALSE,
         color =  colorRampPalette(sequential_hcl(n = 7, palette = "Blues 3"))(64),
         annotation_col = annotation_col,
         annotation_colors =ann_colors
        )})
indata<-scale(LNC[[4]][LNC[[6]],grep("T",colnames(LNC[[4]]))])
F3.3<-as.ggplot(function(){pheatmap(indata[,index],
        cluster_row = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
        fontsize=10,
        legend = FALSE,
        annotation_legend = FALSE,
         color=colorRampPalette(sequential_hcl(n = 7, palette = "Viridis"))(64),
          show_colnames=FALSE ,
         annotation_colors =ann_colors
)})

bubble.df=as.matrix(GetAssayData(COAD1, slot = "counts"))[TOP$gene,]
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
COAD1$CB=rownames(bubble.df)
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,COAD1@meta.data[,c("CB","clusters")],by = "CB")
bubble.df$CB=NULL

celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$clusters)) {
  bubble.df_small=bubble.df%>%filter(clusters==i)
  for (j in TOP$gene) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}

plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
plotdf$celltype=factor(plotdf$celltype,levels = sort(unique(plotdf$celltype)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(TOP$gene)))
if(min(plotdf$exp)<0){plotdf$exp=plotdf$exp+abs(min(plotdf$exp))}
plotdf$exp=ifelse(plotdf$exp>3,3,plotdf$exp)
F3.4<-plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+xlab("Cluster")+
theme_classic()+
  theme(
    axis.text.x.bottom = element_text(hjust = 0, vjust = 0, angle = 0)
  )+theme1
F3.4




setwd("/mnt/RA/whn/crc/our_data/")
QQ<-LNC[[4]][,grep("T",colnames(LNC[[4]]))]%>%scale()
ESTIMATE_score<-read.csv("ESTIMATE_score2.csv",row.names = 1)
ESTIMATE_score<-ESTIMATE_score[paste0("X",colnames(QQ)),]
head(ESTIMATE_score)
#cluster分组各个评分的比较
ESTIMATE_score$cluster = res
ESTIMATE_score$cluster = ifelse(ESTIMATE_score$cluster=="1","cluster1","cluster2")
head(ESTIMATE_score)
# boxplot绘制
my_comparisons <- list( c("cluster1", "cluster2"))
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 12),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position="none",#第一个离Y轴距离，第二十X轴图例在绘图区域的位置
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 10,margin = margin(r=8)),
          #legend.background = element_rect( linetype="solid",colour ="black")
    )
}

# 绘制各种箱型图
p1 <- ggplot(ESTIMATE_score,aes(x=cluster,y=StromalScore),palette = "jco",
             add = "jitter")+xlab("") + ylab("Stromal Score") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$StromalScore)*1.3)

p2 <- ggplot(ESTIMATE_score,aes(x=cluster,y=ImmuneScore),palette = "jco",
             add = "jitter")+xlab("") + ylab("Immune Score") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$ImmuneScore)*1.3)

p3 <- ggplot(ESTIMATE_score,aes(x=cluster,y=ESTIMATEScore),palette = "jco",
             add = "jitter")+xlab("") + ylab("ESTIMATE Score") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$ESTIMATEScore)*1.3)

p4 <- ggplot(ESTIMATE_score,aes(x=cluster,y=TumorPurity),palette = "jco",
             add = "jitter")+xlab("") + ylab("Tumor Purity") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$TumorPurity)*1.3)

p <- p1+p2+p3+p4+plot_layout(nrow = 1)
p
topptx(p,filename = "./plot/ourESTIMATE.pptx",width = 12,height = 3,append = F)




Univariate_cox_fig <- function(td,time,status,
                               boxcol="#b6854d",
                               boxsize=0.15,title = "TCGA"){
    single_cox_fillter <- function(c,time,status){
      cox_data <- as.data.frame(c)
      colnames(cox_data) <- gsub(colnames(cox_data), pattern='-', replacement='_' )  
      pvalue <- c()
      coef <- c()
      HR <- c()
      lower_ci <- c()  
      upper_ci <- c()  
      genes <- colnames(cox_data)[!(colnames(cox_data) %in% c(time, status))]
      for (i in genes) {
        unique_values <- length(table(td[[i]]))
        if (unique_values == 1) {
          pvalue <- c(pvalue, NA)
          coef <- c(coef, NA)
          HR <- c(HR, NA)
          lower_ci <- c(lower_ci, NA)
          upper_ci <- c(upper_ci, NA)
          next
        }
        formula_str <- paste("Surv(", time, ", ", status, ") ~", i)
        s <- as.formula(formula_str)
        model <- coxph(s, data=cox_data)
        sui_cox <- summary(model)
        pvalue <- c(pvalue,sui_cox[["coefficients"]][1,5]) # p值
        coef <- c(coef,cox.zph(model)$table[1,3])   # 回归系数
        # HR <- c(HR,x$conf.int[, 1]) # HR exp(coef)
        HR <- c(HR,exp(sui_cox[["coefficients"]][1, 1]))
        conf_int <- exp(confint(model)) # 置信区间
        lower_ci <- c(lower_ci, conf_int[1, 1])
        upper_ci <- c(upper_ci, conf_int[1, 2])
      }
      cox_count <- data.frame(
        gene_name = genes,
        pvalue = pvalue,
        coef = coef,
        HR = HR,
        lower_ci = lower_ci,
        upper_ci = upper_ci)
      # gene <- cox_count[cox_count$pvalue < 0.005 & cox_count$coef > 0.5,]
      return(cox_count)
    }
    cox_exp <- single_cox_fillter(td,time,status)
    
    ## 单因素cox结果数据整理
    cox_count <- cox_exp
    tabletext <- cbind(cox_count$gene_name,  
                       sprintf("%.2f", cox_count$HR),  
                       sprintf("%.2f", cox_count$lower_ci),
                       sprintf("%.2f", cox_count$upper_ci),  
                       sprintf("%.4f [%.4f, %.4f]",cox_count$HR , cox_count$lower_ci, cox_count$upper_ci),  
                       cox_count$pvalue)  %>% data.frame()
    p <- ifelse(cox_exp$pvalue < 0.001, paste0(sprintf("%.2e", cox_exp$pvalue)," ***"),
                         ifelse(cox_exp$pvalue < 0.01, paste0(sprintf("%.2e", cox_exp$pvalue),"  **"),
                               ifelse(cox_exp$pvalue < 0.05, paste0(round(cox_exp$pvalue, 4),"   *"), paste0(round(cox_exp$pvalue, 4),"    "))))
    tabletext$X6 <- p
    colnames(tabletext) <- c("Characteristics", "exp(coef)", "lower .95", "upper .95", "exp(coef) [confint]", "p")
    tabletext[, 2:4] <- lapply(tabletext[, 2:4], as.numeric)
    print(head(tabletext,2))
    result <- rbind(
      c("Characteristic", NA, NA, NA, "HR (95%CI)", "P"),tabletext)
    result[,2:4] <- lapply(result[,2:4], as.numeric)
    
    ## 画图
    library(forestplot)
    fig <- function(result){
        fig <- forestplot(
            title = title ,  # 添加主标题
          result[, c(1, 5, 6)],   # 选择要在森林图中显示的数据列，第1、5、6列
          mean = result[, 2],     # 指定均值数据列（HR），它将显示为森林图的小方块
          lower = result[, 3],    # 指定95%置信区间的下限数据列
          upper = result[, 4],    # 指定95%置信区间的上限数据列，这些数据将显示为线段穿过方块
          zero = 1,               # 设置零线或参考线为HR=1，这是x轴的垂直线
          fn.ci_norm = "fpDrawDiamondCI",  # 置信区间的形状（这里是钻石形）
          col = fpColors(                  # 颜色设置
            box = boxcol,                 # 方框颜色
            lines = "#8d8680",               # 线条颜色
            zero = "black"                 # 均值为0时的颜色
              
          ),
          boxsize = boxsize,          # 设置小方块的大小
          graphwidth = unit(.34, "npc"),   # 森林图的宽度
          graph.pos = 2 ,          # 指定森林图应该插入到图形中的位置，这里是第2列
          is.summary = c(T, rep(F,nrow(result)-1)),
          txt_gp = fpTxtGp(                # 文本样式的设置
            label = gpar(cex = 1),       # 标签的大小
            ticks = gpar(cex = 1),         # 刻度标记的大小
            xlab = gpar(cex = 1),        # x轴标签的大小
            title = gpar(cex = 1)        # 标题的大小
          )
          # ,line = gpar(lwd = 2) # 修改字体粗细
        )
        return(fig)
    }
    fig <- fig(result)
    return(list(cox_exp,fig))
}
data<-coldata[coldata$origin=="TCGA",]
colnames(data)[3]<-"Cluster"
P1<-Univariate_cox_fig(data[,c("time","vital_status","Cluster")],
                       time="time",status="vital_status" )
					   
					   
#重新绘制热图
library(colorspace) 

#

annCol <- data.frame(
  Cluster = paste0("C",coldata$text),row.names = coldata$match_id)   # 聚成2类的聚类的结果，1，2
annCol$`Tumor Purity` = coldata$TumorPurity

annCol$`Stromal Score` = coldata$StromalScore
annCol$`Immune Score` = coldata$ImmuneScore
annCol$`ESTIMATE Score` = coldata$ESTIMATEScore

annCol$`Tumor Stage`<- as.factor(coldata$tumor_stage)
annCol$Age <- coldata$Age
annCol$Dataset <- coldata$origin
#annCol$OS <- as.factor(data_clin$vital_status)
#annCol$OS.time <- data_clin$time/12
#annCol$Disease <- data2$Primary.Site.Disease
#annCol$Grade <- data2$Tumour.Grade
annColors <- list(Cluster = c("C1" = cols[5],"C2" = cols[6]),
                  `Tumor Purity` =colorRampPalette((c("#f4e7d2","#053061")))(5),                  
                  `Stromal Score` = colorRampPalette((c("#f4e7d2","OrangeRed")))(5),
                  `Immune Score` = colorRampPalette((c("#f4e7d2","MediumVioletRed")))(5),
                  `ESTIMATE Score` = colorRampPalette((c("#f4e7d2","Gold")))(5),
                  Dataset =c("TCGA"=cols[1],"E-MTAB-12862"=cols[2]),
                  `Tumor Stage` = c("1" = "#ffffff","2" = "#88c4e8","3"="#439bd4","4"="#0074b3") ,
                 # Grade = c("N0/1" = "#3CB371","N2/3" = "#9370DB","NX"="#FFA500"),
                 # Disease = c("T1/2" = "#C71585","T3/4" = "#CD5C5C","TX"="#BDB76B"),
                  Age = colorRampPalette((c("#c0cddd","#5b6a76")))(100)
                  
                 # OS =  c("0" = "#EAE3E3", "1" = "#004D61"),
                 # OS.time = colorRampPalette((c("white","#A0522D")))(100))
                  )




#dev.off();
coldata<-coldata[order(coldata$text,coldata$origin),]
library(tidyplots)
plotdata <- as.matrix(dist.r)[coldata$match_id,coldata$match_id]  # 聚成两类的数据
#tiff(filename = "./plot/F3_3.tiff",width = 10,height = 10)
options(repr.plot.width=10, repr.plot.height=10)
F4.1<-pheatmap(mat = plotdata,
        # color = colorRampPalette(sequential_hcl(n = 7, palette = "Blues 3"))(10),
                 #   color=colorRampPalette(rev(colors_discrete_seaside[1:3]))(64),    
         border_color = NA,
         cluster_rows = F,#dendsort(cc[[2]][[2]]),  # 聚成两类的属性，行
         cluster_cols = F,#dendsort(cc[[2]][[2]]),  # 聚成两类的属性，列
         annotation_col = annCol,
         annotation_colors = annColors,
         show_colnames = FALSE,
         show_rownames = FALSE)
		 
		 
		 


## 不同类之间免疫微环境的差异###################
setwd("/mnt/RA/whn/crc/our_data/")
#读取数据，行是基因，列是样本,对保存的文件，要保证我们第一列是基因名，这个需要我们excel修改
in.gct.file = "./ESTIMATE_input.gct"
#####加了GeneSymbol csv转txt
exp.file ="./TCGA_TPM.txt"
#转换成gct文件格式
outputGCT(exp.file, in.gct.file)

filterCommonGenes(input.f = exp.file,
                  output.f = in.gct.file,
                  id = "GeneSymbol")
## start estimate 开始分析
out.score.file = "ESTIMATE_score.gct"
estimateScore(in.gct.file, 
              out.score.file)
#画出每个样本的肿瘤浸润度
plotPurity(out.score.file )
#将结果保存为其他格式
ESTIMATE_score = read.table(out.score.file,
                            skip = 2,
                            header = T,
                            row.names = 1)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[,2:ncol(ESTIMATE_score)]))
ESTIMATE_score$Samples = rownames(ESTIMATE_score)
ESTIMATE_score = ESTIMATE_score[,c(ncol(ESTIMATE_score),2:ncol(ESTIMATE_score)-1)]
#write.xlsx(ESTIMATE_score, "ESTIMATE_score.xlsx",row.names = F)
head(ESTIMATE_score)
write.csv(ESTIMATE_score,"ESTIMATE_score.csv")
setwd("/mnt/RA/whn/crc/our_data/")
index<-order(FSV$text)
QQ<-FSV[index,]
rownames(QQ)<-QQ$match_id
ESTIMATE_score<-read.csv("ESTIMATE_score.csv",row.names = 1)
ESTIMATE_score<-ESTIMATE_score[rownames(QQ),]

#cluster分组各个评分的比较
ESTIMATE_score$cluster = QQ$text
ESTIMATE_score$cluster = ifelse(ESTIMATE_score$cluster=="1","cluster 1","cluster 2")
head(ESTIMATE_score)
# boxplot绘制
my_comparisons <- list( c("cluster 1", "cluster 2"))
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),#X轴标签加粗及字体大小     
          axis.text.y = element_text(face = "bold",size = 12),#y坐标轴刻度标签加粗
          axis.text.x = element_text(face = "bold",size = 12),#x坐标轴刻度标签加粗 
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.position="none",#第一个离Y轴距离，第二十X轴图例在绘图区域的位置
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 10,margin = margin(r=8)),
          #legend.background = element_rect( linetype="solid",colour ="black")
    )
}

# 绘制各种箱型图
p1 <- ggplot(ESTIMATE_score,aes(x=cluster,y=StromalScore),palette = "jco",
             add = "jitter")+xlab("") + ylab("Stromal Score") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$StromalScore)*1.3)

p2 <- ggplot(ESTIMATE_score,aes(x=cluster,y=ImmuneScore),palette = "jco",
             add = "jitter")+xlab("") + ylab("Immune Score") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$ImmuneScore)*1.3)

p3 <- ggplot(ESTIMATE_score,aes(x=cluster,y=ESTIMATEScore),palette = "jco",
             add = "jitter")+xlab("") + ylab("ESTIMATE Score") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$ESTIMATEScore)*1.3)

p4 <- ggplot(ESTIMATE_score,aes(x=cluster,y=TumorPurity),palette = "jco",
             add = "jitter")+xlab("") + ylab("Tumor Purity") +
  geom_boxplot(aes(fill=cluster),position=position_dodge(0.5),width=0.6)+
  scale_fill_manual(values = c(cols[5],cols[6])) +
  theme_zg() +
  stat_compare_means(comparisons = my_comparisons,method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001,
 0.01, 0.05, Inf),symbols = c("****", "***", "**", "*",  "ns"))) +
  stat_compare_means(label.y = max(ESTIMATE_score$TumorPurity)*1.3)

p <- p1+p2+p3+p4+plot_layout(nrow = 1)
p
topptx(p,filename = "./plot/ESTIMATE.pptx",width = 12,height = 3,append = F)