#清理内存
#rm(list = ls())
#设置路径
#setwd("D:/深圳市人民医院/PAGE")
#加载包
library(gtools)
library(stringr)
#install.packages("snowfall")
library(snowfall)
### 读取表达谱 清空无关变量
# EXP<-COAD_methylation
# Clinic<-group_inf
# cutoff=0.001
# C=3

GeneID_2_Symbol<-function(Result){
  COL1<-as.matrix(coding_gene[match(Result[,1],coding_gene$entrez_id),]$symbol)
  COL2<-as.matrix(coding_gene[match(Result[,2],coding_gene$entrez_id),]$symbol)
  TABLE<-cbind(COL1,Result,COL2)
  colnames(TABLE)<-c("Gene Symbol","Gene ID","Gene ID","Gene Symbol")
  return(TABLE)
}
REOF1<-function(EXP,Clinic,n,cutoff,C){
  ####分类###############
  Clinic<-Clinic
  
  ######## 获取差异基因######################
  
  Result<-GetDEG(EXP,Clinic,n)
  # Result[[1]]<-DEG[[1]]    A中的差异基因
  # Result[[2]]<-DEG[[2]]    B中的差异基因
  # Result[[3]]<-DEG[[3]]    AB 共有的差异基因
  # Result[[4]]<-EXP_TypeA   差异前的A表达谱
  # Result[[5]]<-EXP_TypeB   差异前的B表达谱
  
  ############# 获取差异表达谱###############
  
  Mat1<-GetDEG_exp(Result[[3]],Result[[4]],RowN = T)
  Mat2<-GetDEG_exp(Result[[3]],Result[[5]],RowN = T)
  
  rownames(Result[[5]])
  ############稳定对筛选方法################ 
  ############## Fisher
  Stable_pair_fisher<-function(F){
    GeneID<-rownames(Mat1)
    Nsample<-ncol(Mat1)
    Result<-matrix(NA,ncol=3)
    pb <- txtProgressBar(style=3)
    for( i in 1:length(F)){
      i<-F[i]
      CopyMat<-matrix(rep(Mat1[i,],times=(nrow(Mat1)-i)),ncol=ncol(Mat1),byrow=T)   #扩充矩阵
      IF_Positive<-CopyMat-Mat1[-(1:i),]  #A:第i行的基因
      # print(IF_Positive)
      IF_Positive[which(IF_Positive<0)]<- 0   #
      IF_Positive[which(IF_Positive>0)]<- 1
      Sum_AmB<-rowSums(IF_Positive)             #表达值A大于B的样本个数
      temp<-as.matrix(rep(dim(IF_Positive)[2],times = dim(IF_Positive)[1]))-Sum_AmB
      Sum_AmB<-cbind(Sum_AmB,temp)
      
      CopyMat<-matrix(rep(Mat2[i,],times=(nrow(Mat2)-i)),ncol=ncol(Mat2),byrow=T)   #扩充矩阵
      IF_Positive<-CopyMat-Mat2[-(1:i),]  #A:第i行的基因
      # print(IF_Positive)
      IF_Positive[which(IF_Positive<0)]<- 0   #
      IF_Positive[which(IF_Positive>0)]<- 1
      Sum_AmB2<-rowSums(IF_Positive)             #表达值A大于B的样本个数
      temp<-as.matrix(rep(dim(IF_Positive)[2],times = dim(IF_Positive)[1]))-Sum_AmB2
      Sum_AmB2<-cbind(Sum_AmB2,temp)
      
      SUM_AMB<-cbind(Sum_AmB,Sum_AmB2)
      
      P<-apply(SUM_AMB,1,function(ind){
        Temp<- matrix(c(ind[1], ind[3],ind[2],ind[4]), nrow = 2)
        Temp2<- fisher.test(Temp)
        return(Temp2[1])} )
      P<-as.matrix(p.adjust(as.matrix(unlist(P)),'BH'))
      Gene<-GeneID[-(1:i)]
      TempB<-as.matrix(Gene[which(P<cutoff)])
      if(length(TempB)!=0){
        TempA<-as.matrix(rep(GeneID[i],times=(length(TempB))))
        R_Mat_p1<-cbind(cbind(TempA,TempB),P[which(P<cutoff)])  #B<A
        TempA2<-matrix(SUM_AMB[which(P<cutoff),],ncol = 4)
        
        if(sum(TempA2[,1]>TempA2[,2])!=0)
        {R_Mat_p1[which(TempA2[,1]>TempA2[,2]),1:2]<-R_Mat_p1[which(TempA2[,1]>TempA2[,2]),2:1]}
      }
      else{
        R_Mat_p1<-matrix(NA,ncol=3)
      }      
      Result<-rbind(Result,R_Mat_p1)
      setTxtProgressBar(pb, i/nrow(Mat1))
    } 
    return(na.omit(Result))   #A>B
  }
  Term<-function(n,c){
    
    A<-matrix(c(1:(n-1)))
    B<-rep(c(1:c,c:1),times=(n/c+1))
    rownames(A)<-B[c(1:(n-1))]
    C<-list()
    for(i in 1:c){
      C[[i]]<-c(A[which(rownames(A)==i),])
    }
    return(C)
  }
  
  #BiocManager::install("snowfall")
  library(snowfall)
  sfStop()
  sfInit(parallel = TRUE, cpus = C)
  sfExport("Mat1", "Mat2","cutoff") 
  result <- sfLapply(Term(length(Result[[3]]),C), Stable_pair_fisher)
  sfStop()
  RES<-do.call("rbind",result)
  #############   阈值  
  
  # R_A<-Stable_pair2(DEG_EXP_A,cutoff)           #Stable_pair2
  # R_B<-Stable_pair2(DEG_EXP_B,cutoff)
  
  #######################################################    
  
  print("The number of Reverse_Pair")
  print(length(RES)/3)
  RES<-RES[order(as.numeric(RES[,3])),]
  ##############返回需要的信息#################
  WHN<-list()
  WHN[[1]]<-RES
  return(WHN)  
}
Get_Clinic_label<-function(Clinic){
  Temp<-as.matrix(rep(NA,times=(dim(Clinic)[1])))
  #####自定义####
  Temp[which(Clinic$condition==1)]<-1
  Temp[which(Clinic$condition==0)]<-2
  colnames(Temp)<-"label"
  Clinic<-cbind(Clinic,Temp)
  return(Clinic)
} 
GetDEG<-function(EXP,Clinic,n){  #clinic 统一通过lable中的0和1判断   输入的EXP是有行名
  
  TypeA<-Clinic[which(Clinic$label==1),1]#回满足条件Clinic$label==1的行的索引,并且只选择第一列的数据
  TypeB<-Clinic[which(Clinic$label==2),1]
  print("Sample of TypeA")         #### A类名称
  print(length(TypeA))
  print("Sample of TypeB")         #### B类名称
  print(length(TypeB))
  EXP<-EXP
  
  EXP_TypeA<-EXP[,match(TypeA,colnames(EXP))]    ####  1对应于原GBM  0对应与LGG
  EXP_TypeB<-EXP[,match(TypeB,colnames(EXP))]
  
  Rank_TypeB<-GetRank(EXP_TypeB)
  Rank_TypeA<-GetRank(EXP_TypeA)
  
  
  MedRank_B<-MedRank(Rank_TypeB)
  MedRank_A<-MedRank(Rank_TypeA)
  
  DEG<-SelectDEG(MedRank_A,MedRank_B,n)
  Result<-list()
  Result[[1]]<-DEG[[1]]
  Result[[2]]<-DEG[[2]]
  Result[[3]]<-DEG[[3]]
  Result[[4]]<-EXP_TypeA
  Result[[5]]<-EXP_TypeB
  
  return(Result)
}
GetRowNames<-function(EXP){
  
  
  NAMES<-EXP[,1]
  Cnames<-colnames(EXP)
  EXP<-EXP[,-1]
  rownames(EXP)<-NAMES
  colnames(EXP)<-Cnames[-1]
  return(as.matrix(EXP))
}
#####基因ID和numb#######
GeneID_2_Symbol<-function(Result){
  COL1<-as.matrix(coding_gene[match(Result[,1],coding_gene$entrez_id),]$symbol)
  COL2<-as.matrix(coding_gene[match(Result[,2],coding_gene$entrez_id),]$symbol)
  TABLE<-cbind(COL1,Result,COL2)
  colnames(TABLE)<-c("Gene Symbol","Gene ID","Gene ID","Gene Symbol")
  return(TABLE)
}
#获取秩次矩阵
#EXP_rank<-EXP_TypeB
GetRank<-function(EXP_rank){ 
  
  N<-ncol(EXP_rank)
  NR<-nrow(EXP_rank)
  for(i in (1:NR)){
    Temp<-length(which(EXP_rank[i,]==0))     ##计算0值数量
    if(Temp>(N/2)){
      EXP_rank[i,]<-NA
    }
  }
  EXP_rank<-na.omit(EXP_rank)           #删除0值多的基因
  NewMat<-rownames(EXP_rank)
  #print(head(NewMat))###基因ID
  for (i in (1:N)){
    Temp1<-unique(as.numeric(EXP_rank[,i]))    #样本所以基因表达值 去重
    numb<-length(Temp1)                        #去重后还剩多少基因
    OrderSqe<-matrix(data=(1:numb),ncol=1)     #生成序列数列
    EXP_Drd<-cbind(as.matrix(sort(Temp1)),OrderSqe)    #排序
    Rank<-EXP_Drd[match(EXP_rank[,i],EXP_Drd[,1]),2]   #匹配
    NewMat<-cbind(NewMat,Rank)             #得到秩序[]
    
  }
  NewMat<-GetRowNames(NewMat)
  
  colnames(NewMat)<-colnames(EXP_rank)
  return(NewMat)   #返回表达谱的基因秩序  S*N
}

##取中位秩次
MedRank<-function(RankMat){
  Temp<-as.matrix(apply(RankMat,1,function(x){median(as.numeric(x))}))
  return(Temp)
}
# RankA<-MedRank_A
# RankB<-MedRank_B
SelectDEG<-function(RankA,RankB,n){      ##n为选取秩序差异最大的2n个基因   rank n*2

  GeneID_A<-as.matrix(rownames(RankA))
  GeneID_B<-as.matrix(rownames(RankB))
  RankA<-cbind(GeneID_A,RankA)
  RankB<-cbind(GeneID_B,RankB)           #将行名改为第一列
  #print("n of RA")
  #print(length(RankA))
  #print(head(RankA))
  #print("n of RB")
  #print(length(RankB))
  #print(head(RankB))
  
  InterGeneID<-intersect(GeneID_A,GeneID_B)##交叠的基因
  #print("ING")
  #print(head(InterGeneID))
  RankA<-RankA[match(InterGeneID,RankA[,1]),]
  RankB<-RankB[match(InterGeneID,RankB[,1]),]          ##取交叠基因的秩序矩阵
  
  Def_Rank<-as.matrix(as.numeric(RankA[,2])-as.numeric(RankB[,2]))             ##秩序差   ####   eg.LGG-GBM   疾病-正常
  rownames(Def_Rank)<-InterGeneID
  
  A<-names(Def_Rank[order(as.numeric(Def_Rank),decreasing = F),])[c(1:n)] 
  C<-names(Def_Rank[order(as.numeric(Def_Rank),decreasing = T),])[c(1:n)] 
  
  DEG<-list()
  if(length(intersect(A,C))==0){
    DEG[[1]]<-A
    DEG[[2]]<-C
  }else{
    Overlop<-intersect(A,C)
    DEG[[1]]<-setdiff(A,Overlop)      #去掉交叠的基因
    DEG[[2]]<-setdiff(C,Overlop)
  }
  
  DEG[[3]]<-unique(as.matrix(c(A,C)))
  print(dim(DEG[[3]]))
  names(DEG)<-c("Gene1","Gene2","All_DEG")
  #print("n of All_DEG")
  
  return(DEG)
  
}
## 获取差异基因表达谱
GetDEG_exp<-function(DEG,EXP,RowN = T)
{
  if(RowN==T){
    GeneID<-as.matrix(rownames(EXP))
  }else{
    GeneID<-as.matrix(EXP[,1])
  }
  
  EXP[match(DEG,GeneID),] 
  
}

ScoreMat2<-function(DEG_pair,EXP,N
=-1){   #DEG_infor为训练出的差异基因信息矩阵

DEG_pair<-DEG_pair[[1]]
if(length(DEG_pair)==3){
DEG_pair<-as.data.frame(t(DEG_pair))
}
DEG<-unique(c(DEG_pair[,1],DEG_pair[,2]))

sub_EXP<-GetDEG_exp(DEG,EXP,RowN = T)
scoreMat<-matrix(NA,nrow=nrow(DEG_pair),ncol=ncol(sub_EXP))

for (i in 1:nrow(DEG_pair)){

        EXP_A_B<-sub_EXP[DEG_pair[i,c(1:2)],]
        Temp<-t(as.matrix(EXP_A_B[1,]-EXP_A_B[2,]))   #相减  A-B
        Temp[which(Temp>0)]<-1     #A>B
        Temp[which(Temp<=0)]<-N   #A<B
        scoreMat[i,]<-Temp
    }
   
    colnames(scoreMat)<-colnames(sub_EXP)
    rownames(scoreMat)<-c(1:nrow(scoreMat))
    return(na.omit(scoreMat))
}



