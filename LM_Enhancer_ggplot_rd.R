setwd("/Users/mar/Research/01_PhD/mm10/13_PredictiveModels/06_PredictiveModelsDiff/ROUND3")

library(caret)
library(viridis)
library(ggpubr)

# Load and prepare data
FileMES<-read.table("gene_PE_H3K4me1_H3K4me3_H3K36me3_H3K27me3_H3K27ac_expression_FDR0_P100_MES.txt")
dataMES<-cbind(FileMES[,3:8])
colnames(dataMES)<-c("H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K27ac","expression")
dataMES<-log2(dataMES+0.1)
dataMES<-data.frame(dataMES)

FileCP<-read.table("gene_PE_H3K4me1_H3K4me3_H3K36me3_H3K27me3_H3K27ac_expression_FDR0_P100_CP.txt")
dataCP<-cbind(FileCP[,3:8])
colnames(dataCP)<-c("H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K27ac","expression")
dataCP<-log2(dataCP+0.1)
dataCP<-data.frame(dataCP)

FileCM<-read.table("gene_PE_H3K4me1_H3K4me3_H3K36me3_H3K27me3_H3K27ac_expression_FDR0_P100_CM.txt")
dataCM<-cbind(FileCM[,3:8])
colnames(dataCM)<-c("H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K27ac","expression")
dataCM<-log2(dataCM+0.1)
dataCM<-data.frame(dataCM)

FileNPC<-read.table("gene_PE_H3K4me1_H3K4me3_H3K36me3_H3K27me3_H3K27ac_expression_FDR0_P100_NPC.txt")
dataNPC<-cbind(FileNPC[,3:8])
colnames(dataNPC)<-c("H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K27ac","expression")
dataNPC<-log2(dataNPC+0.1)
dataNPC<-data.frame(dataNPC)

FileCN<-read.table("gene_PE_H3K4me1_H3K4me3_H3K36me3_H3K27me3_H3K27ac_expression_FDR0_P100_CN.txt")
dataCN<-cbind(FileCN[,3:8])
colnames(dataCN)<-c("H3K4me1","H3K4me3","H3K36me3","H3K27me3","H3K27ac","expression")
dataCN<-log2(dataCN+0.1)
dataCN<-data.frame(dataCN)

# Create random dataTrain
set.seed(24)
expression_rd <- sample(dataMES$expression)
dataMES_rd <- dataMES
dataMES_rd$expression <- expression_rd

expression_rd <- sample(dataCP$expression)
dataCP_rd <- dataCP
dataCP_rd$expression <- expression_rd

expression_rd <- sample(dataCM$expression)
dataCM_rd <- dataCM
dataCM_rd$expression <- expression_rd

expression_rd <- sample(dataNPC$expression)
dataNPC_rd <- dataNPC
dataNPC_rd$expression <- expression_rd

expression_rd <- sample(dataCN$expression)
dataCN_rd <- dataCN
dataCN_rd$expression <- expression_rd

# 10X cross-validation

# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats = 3)

# learn model

modelMES <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                  data=dataMES, trControl=train_control, method="lm")
modelCP <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                 data=dataCP, trControl=train_control, method="lm")
modelCM <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                 data=dataCM, trControl=train_control, method="lm")
modelNPC <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                  data=dataNPC, trControl=train_control, method="lm")
modelCN <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                 data=dataCN, trControl=train_control, method="lm")

modelMES_rd <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                     data=dataMES_rd, trControl=train_control, method="lm")
modelCP_rd <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                    data=dataCP_rd, trControl=train_control, method="lm")
modelCM_rd <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                    data=dataCM_rd, trControl=train_control, method="lm")
modelNPC_rd <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                     data=dataNPC_rd, trControl=train_control, method="lm")
modelCN_rd <- train(expression~H3K4me1+H3K4me3+H3K36me3+H3K27me3+H3K27ac,
                    data=dataCN_rd, trControl=train_control, method="lm")

# define funtion to visualize model performance/results
results <- function(model,name,dataMES,dataCP,dataCM,dataNPC,dataCN,model_name){
  # summarize results
  model_sum<-capture.output(summary(model))
  out_sum = paste("outputs/",name,"_summary.txt", sep = "")
  writeLines(model_sum,out_sum)
  model_print<-capture.output(print(model))
  out_print = paste("outputs/",name,"_print.txt", sep = "")
  writeLines(model_print,out_print)
  model_imp<-capture.output(varImp(model))
  out_imp = paste("outputs/",name,"_varImp.txt", sep = "")
  writeLines(model_imp,out_imp)
  
  out_varImp<-paste("plots/",name,"_varImp.pdf",sep="")
  p<-ggplot(varImp(model,scale = FALSE))+
    geom_bar(stat="identity", fill="darkviolet")+
    labs(title="Enhancer model variable importance",x="", y = "importance") +
    theme_bw() +
    theme(legend.position="bottom",axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"),
          plot.title = element_text(size=15,face="bold"),strip.text = element_text(size=12,face="bold")) 
  ggsave(out_varImp,plot=p,device = "pdf",width = 5, height = 4)
  
  # plot predicted vs. measured
  exp_predMES<-predict(model,dataMES)
  exp_predCP<-predict(model,dataCP)
  exp_predCM<-predict(model,dataCM)
  exp_predNPC<-predict(model,dataNPC)
  exp_predCN<-predict(model,dataCN)
  #model_res <- lm(exp_pred~expression, data=dataMES)
  #summary(model_res)
  r_MES<-cor(exp_predMES,dataMES$expression)
  r_CP<-cor(exp_predCP,dataCP$expression)
  r_CM<-cor(exp_predCM,dataCM$expression)
  r_NPC<-cor(exp_predNPC,dataNPC$expression)
  r_CN<-cor(exp_predCN,dataCN$expression)
  
  predicted<-c(exp_predMES,exp_predCP,exp_predCM,exp_predNPC,exp_predCN)
  measured<-c(dataMES$expression,dataCP$expression,dataCM$expression,dataNPC$expression,dataCN$expression)
  Cell<-c(rep("MES",nrow(dataMES)),rep("CP",nrow(dataCP)),rep("CM",nrow(dataCM)),rep("NPC",nrow(dataNPC)),rep("CN",nrow(dataCN)))
  c<-data.frame(predicted,measured,Cell)
  c$Cell = factor(c$Cell, levels = c("MES","CP","CM","NPC","CN"))
  
  text<-data.frame(label=c(paste("r = ",round(r_MES,2)),paste("r = ",round(r_CP,2)),paste("r = ",round(r_CM,2)),
                           paste("r = ",round(r_NPC,2)),paste("r = ",round(r_CN,2))),
                   Cell=c("MES","CP","CM","NPC","CN"))
  out_plots<-paste("plots/",name,".pdf",sep="")
  p<-ggplot(c) +
    geom_hex(aes(measured, predicted), bins = 100) +
    scale_fill_gradientn("", colours = rev(viridis(300)))+
    geom_smooth(aes(measured, predicted),method = "lm",level=0)+
    labs(title=paste("Expression prediction using",model_name),x="Measured expression (log(FPKMs + 0.1))", y = "Predicted expression") +
    geom_text(data = text, mapping = aes(x = -Inf, y = -Inf, label = label),hjust = -2, vjust = -1, size=7) +
    theme_bw() +
    theme(legend.position="right",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
          strip.text = element_text(size=20,face="bold")) +
    facet_wrap(~ Cell, scales="free")
  ggsave(out_plots,plot=p,device = "pdf",width = 15, height = 10)
  
  r<-c(r_MES,r_CP,r_CM,r_NPC,r_CN)
  r
}


r_MES<-results(modelMES,"LM_E_MES",dataMES,dataCP,dataCM,dataNPC,dataCN,"MES model")
r_CP<-results(modelCP,"LM_E_CP",dataMES,dataCP,dataCM,dataNPC,dataCN,"CP model")
r_CM<-results(modelCM,"LM_E_CM",dataMES,dataCP,dataCM,dataNPC,dataCN,"CM model")
r_NPC<-results(modelNPC,"LM_E_NPC",dataMES,dataCP,dataCM,dataNPC,dataCN,"NPC model")
r_CN<-results(modelCN,"LM_E_CN",dataMES,dataCP,dataCM,dataNPC,dataCN,"CN model")

r_MES_rd<-results(modelMES_rd,"LM_E_MES_rd",dataMES,dataCP,dataCM,dataNPC,dataCN,"MES model")
r_CP_rd<-results(modelCP_rd,"LM_E_CP_rd",dataMES,dataCP,dataCM,dataNPC,dataCN,"CP model")
r_CM_rd<-results(modelCM_rd,"LM_E_CM_rd",dataMES,dataCP,dataCM,dataNPC,dataCN,"CM model")
r_NPC_rd<-results(modelNPC_rd,"LM_E_NPC_rd",dataMES,dataCP,dataCM,dataNPC,dataCN,"NPC model")
r_CN_rd<-results(modelCN_rd,"LM_E_CN_rd",dataMES,dataCP,dataCM,dataNPC,dataCN,"CN model")

c<-data.frame(r = c(r_MES[-1],r_CP[-2],r_CM[-3],r_NPC[-4],r_CN[-5],
                    r_MES_rd[-1],r_CP_rd[-2],r_CM_rd[-3],r_NPC_rd[-4],r_CN_rd[-5]),
              Cell = c(rep("MES",4),rep("CP",4),rep("CM",4),rep("NPC",4),rep("CN",4),
                       rep("MES",4),rep("CP",4),rep("CM",4),rep("NPC",4),rep("CN",4)),
              Method = c(rep("Model",20),rep("Random\nmodel",20)))
c$Cell = factor(c$Cell, levels = c("MES","CP","CM","NPC","CN"))

p<-ggplot(c, aes(x=Method, y=r,fill=Method)) + 
  scale_fill_manual(values=c("cornflowerblue","seagreen")) +
  geom_boxplot()+
  labs(title="Comparison between the performances of the models and the random models",x="", y = "r") +
  theme_bw() +
  theme(legend.position="bottom",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  stat_compare_means(paired = TRUE, method = "t.test",aes(label = ..p.signif..),label.x = 1.5,vjust = 1, size = 10) +
  facet_wrap(~ Cell, scales="free")
ggsave("plots/r_E_boxplot.pdf",plot=p,device = "pdf",width = 15, height = 10)





shapiro.test(r_MES[-1])
shapiro.test(r_MES_rd[-1])
shapiro.test(r_CP[-2])
shapiro.test(r_CP_rd[-2])
shapiro.test(r_CM_rd[-3])
shapiro.test(r_CM[-3])
shapiro.test(r_NPC[-4])
shapiro.test(r_NPC_rd[-4])
shapiro.test(r_CN_rd[-5])
shapiro.test(r_CN[-5])
shapiro.test(c(r_MES_rd[-1],r_MES[-1]))
shapiro.test(c(r_CP_rd[-2],r_CP[-2]))
shapiro.test(c(r_CM_rd[-3],r_CM[-3]))
shapiro.test(c(r_NPC_rd[-4],r_NPC[-4]))
shapiro.test(c(r_CN_rd[-5],r_CN[-5]))

