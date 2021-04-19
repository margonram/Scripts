setwd("/Users/mgonzalez/Research/01_PhD/mm10/13_PredictiveModels/01_PredictiveModelsHiC-all/ROUND3")

library(caret)
library(viridis)

# Load and prepare data
File<-read.table("gene_E-H4K20me3-H2Bub1-H3K79me2-H3K36me3-H3K27me2-H3K27me1-H3K4me3-H3K4me1-H3K27ac-H3K27me3_expression_HiC-all.txt")
data<-cbind(File[,3:13])
colnames(data)<-c("H4K20me3","H2Bub1","H3K79me2","H3K36me3","H3K27me2","H3K27me1","H3K4me3","H3K4me1","H3K27ac","H3K27me3","expression")
data<-log2(data+0.1)
data<-data.frame(data)
data$gene<-File[,1]

# Create sets of training and test
set.seed(24)
trainIndex <- createDataPartition(data$expression, p = .8,list = FALSE,times = 1)
dataTrain <- data[ trainIndex,]
dataTest  <- data[-trainIndex,]

# Create random dataTrain
set.seed(11)
expression_rd <- sample(dataTrain$expression)
dataTrain_rd <- dataTrain
dataTrain_rd$expression <- expression_rd

# 10X cross-validation

# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats = 3)

# learn model

model <- train(expression~H4K20me3+H2Bub1+H3K79me2+H3K36me3+H3K27me2+H3K27me1+H3K4me3+H3K4me1+H3K27me3+H3K27ac,
                data=dataTrain, trControl=train_control, method="lm")

model_rd <- train(expression~H4K20me3+H2Bub1+H3K79me2+H3K36me3+H3K27me2+H3K27me1+H3K4me3+H3K4me1+H3K27me3+H3K27ac,
                data=dataTrain_rd, trControl=train_control, method="lm")

# summarize results
model_sum<-capture.output(summary(model))
out_sum = paste("outputs/Enhancer_summary.txt", sep = "")
writeLines(model_sum,out_sum)
model_print<-capture.output(print(model))
out_print = paste("outputs/Enhancer_print.txt", sep = "")
writeLines(model_print,out_print)
model_imp<-capture.output(varImp(model))
out_imp = paste("outputs/Enhancer_varImp.txt", sep = "")
writeLines(model_imp,out_imp)

out_varImp<-paste("plots/Enhancer_varImp.pdf",sep="")
p<-ggplot(varImp(model,scale = FALSE))+
  geom_bar(stat="identity", fill="darkviolet")+
  labs(title="Enhancer model variable importance",x="", y = "importance") +
  theme_bw() +
  theme(legend.position="bottom",axis.text=element_text(size=12,face="bold"),axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size=15,face="bold"),strip.text = element_text(size=12,face="bold")) 
ggsave(out_varImp,plot=p,device = "pdf",width = 5, height = 4)


# plot predicted vs. measured
exp_pred<-predict(model,dataTest)
exp_pred_rd<-predict(model_rd,dataTest)
r<-round(cor(exp_pred,dataTest$expression),2)
r_rd<-round(cor(exp_pred_rd,dataTest$expression),2)

predicted<-c(exp_pred,exp_pred_rd)
measured<-c(dataTest$expression,dataTest$expression)
models<-c(rep("Model",nrow(dataTest)),rep("Random model",nrow(dataTest)))

text<-data.frame(label=c(paste("r = ",r),paste("r = ",r_rd)),models=c("Model","Random model"))
c<-data.frame(predicted,measured,models)

pdf("plots/Enhancer_prediction.pdf", width = 10,height = 5)
ggplot(c) +
  geom_hex(aes(measured, predicted), bins = 100) +
  scale_fill_gradientn("", colours = rev(viridis(300)))+
  geom_smooth(aes(measured, predicted),method = "lm",level=0)+
  labs(title="Expression prediction in the test subset",x="Measured expression (log(FPKMs + 0.1))", y = "Predicted expression") +
  geom_text(data = text, mapping = aes(x = -Inf, y = -Inf, label = label),hjust = -2, vjust = -1, size=7) +
  theme_bw() +
  theme(legend.position="right",axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),legend.title=element_text(size=20),plot.title = element_text(size=24,face="bold"),
        strip.text = element_text(size=20,face="bold")) +
  facet_wrap(~ models)
dev.off()
