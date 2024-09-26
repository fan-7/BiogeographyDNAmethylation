## 20191119
## GAB1000 MEthylation geography inference
## extract sample1000 snp(hwp0.001 exclude) to get PC5 
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/ancestry/GAB_1000G/2509GAB1000G/GAB2509Merge1000G/GAB2509.exHWP --keep /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/GABmale990.samplelist.txt --make-bed --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/GABmale990Methy.genome
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/GAB1000Methy.genome --pca 5 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/pc5.GAB1000genome
 

plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/GABmale990Methy.genome --pca 10 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/pc.GABmalegenome

cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS 
##20191128 manage data for logit
R
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS")
male990 <- read.table("GABmale990Methy.txt", header = TRUE)
pc10 <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/GAB1000genome/pc10.GABmalegenome.eigenvec", header = FALSE)
names(pc10) <- c("FID", "IID", "PC1", "PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
dataGeo <- merge(male990, pc10[,-1], by = "IID")
library(dplyr)
dataGeo <- select(dataGeo, FID, everything())
write.table(na.omit(dataGeo), file = "dataGeo.GAB922", quote = FALSE, row.names = FALSE)
dataGeo <- na.omit(dataGeo)
save(dataGeo, file = "dataGeo.RData")
## use preparer data to EWAS
#
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS")
load("dataGeo.RData")
load("dataMethy.RData")
gc()
blogit <- function(x)
{
# logistic model
fit <- glm(dataGeo$y1 ~ x + dataGeo$age+dataGeo$PC1+dataGeo$PC2+dataGeo$PC3+dataGeo$PC4+dataGeo$PC5, family=binomial())
s <- summary(fit)
blogitRe <- coef(s)[2, c(1,2,3,4)]
return(blogitRe)
}
##error from x length differ from dataGeo
 which(colnames(usedMethy) %in% dataGeo$IID)
 dataMethy <- usedMethy[,which(colnames(usedMethy) %in% dataGeo$IID)]
 ## keep dataGeo and dataMethy dim equle
y1_dataMethy_bligitResult <- apply(dataMethy, 1, blogit)
write.table(y1_dataMethy_bligitResult, file = "y1_dataMethy_bligitResult.Result", quote = F, row.names = F)
## try to prepare code for conditional EWAS 
##!!! but test y1 ~ cpg +age + pc1 -pc10  get P value mean = 0.99
for(i in 1: 20){
    if(i ==1){
	   #filter P < 1e-6 save 
	   tmp <- paste0("awk '$4<1e-6 {print $2,$1,$3,$12}' ",outpath,outname,"_",i,".assoc.logistic"," > ",outpath,resultname,"_",i)
	   system(tmp)
	   # plotManhattan resultfile= resultname_i
	   resultfile <- paste0(resultname,"_",i)
	   dataManhattan <- fread(resultfile, header = FALSE)
	   names(dataManhattan) <- c("SNP", "CHR", "BP", "P")
	   png(paste0(resultfile,".manhattan.png"))
	   manhattan(dataManhattan, col = c("blue4", "red3"),  main = paste0(resultfile,".ManhattanPlot"))
       dev.off()
	   # select top snpfile = resultname_i.selectSNP_i(snplist).
	   selectSNP[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(dataManhattan$P)),])
       rm(dataManhattan)
       gc()
	   covarname <- as.character(selectSNP[i,1])
	   snpfile <- paste0(resultfile,".selectSNP_",i)
	   write.table(covarname, file = snpfile, quote = FALSE, row.names = FALSE)
	   # recode top SNP to raw 012file
	   tmp <- paste0("plink",genofile," --extract ",outpath,snpfile," --recodeA"," --out ",outpath,outname,"_",i)
	   system(tmp)
	   # covar file 
	   covarfile <- read.table(paste0(outname,"_",i,".raw"),header = TRUE)
	   covarname <- names(covarfile)[-c(1:6)]
	}else{
	   # plink covar(top snp)
	   tmp <- paste0("plink",genofile,phenofile,phename," --covar ",outpath,outname,"_",i-1,".raw"," --covar-name ",covarname," --logistic hide-covar --ci 0.95"," --out ",outpath,outname,"_",i)
       system(tmp)
	   #filter P < 1e-6 save 
	   tmp <- paste0("awk '$12<1e-6 {print $2,$1,$3,$12}' ",outpath,outname,"_",i,".assoc.logistic"," > ",outpath,resultname,"_",i)
	   system(tmp)
	   # plotManhattan resultfile= resultname_i
	   resultfile <- paste0(resultname,"_",i)
	   dataManhattan <- fread(resultfile, header = FALSE)
	   # if no snp P <1e-6 break
	   if(nrow(dataManhattan)==0){
	               break
		               }
	   names(dataManhattan) <- c("SNP", "CHR", "BP", "P")
	   png(paste0(resultfile,".manhattan.png"))
	   manhattan(dataManhattan, col = c("blue4", "red3"),  main = paste0(resultfile,".ManhattanPlot"))
       dev.off()
	   # select top snpfile = resultname_i.selectSNP_i(snplist).
	   selectSNP[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(dataManhattan$P)),][1,])
       rm(dataManhattan)
       gc()
	   covarname <- as.character(selectSNP[1:i,1])
	   snpfile <- paste0(resultfile,".selectSNP_",i)
	   write.table(covarname, file = snpfile, quote = FALSE, row.names = FALSE)
	   # recode top SNP to raw 012file
	   tmp <- paste0("plink",genofile," --extract ",outpath,snpfile," --recodeA"," --out ",outpath,outname,"_",i)
	   system(tmp)
	   # covar file 
	   covarfile <- read.table(paste0(outname,"_",i,".raw"),header = TRUE)
	   covarname <- names(covarfile)[-c(1:6)]
	    # replace X to " ", . to ":"(problem from read.table 8: -> X8.)
	   covarname <- gsub("X", " ", covarname)
	   covarname <- gsub("\\.", ":", covarname)
	   covarname <- paste0(covarname, collapse = ",") 
	}
write.table(selectSNP, file = paste0(outpath,"selectSNP.",outname), quote = FALSE, row.names = FALSE)
}
## prepare to plot ManhattanPlot but fail 
library(qqman)
png("yHan_usedMethy_age_bligitResult.qqplot.png")
qq(yHan_usedMethy_age_bligitResult[3, ], main = "y1_usedMethy_bligitResult p-values")
dev.off()
load("usedAno.RData")
gc()
methyID_logitResult <- data.frame(colnames(yHan_usedMethy_age_bligitResult),t(yHan_usedMethy_age_bligitResult))
usedManhattan <- merge(usedAno, methyID_logitResult, by.x ="probID", by.y = "colnames.yHan_usedMethy_age_bligitResult.")
usedManhattanplot <- data.frame(usedManhattan$probID, as.integer(usedManhattan$CHR), usedManhattan$BP, usedManhattan$Pr...z..)
names(usedManhattanplot) <- c("cpg", "CHR", "BP", "P")
gc()
png("yHan_usedMethy_age_bligitResult.manhattanplot.png")
manhattan(usedManhattanplot, col = c("blue4", "red3"), cex = 0.6, cex.axis = 0.9, ymax = 12, main = "yHan_usedMethy_age_bligitResult.manhattanplot")
dev.off()
#####################################################################
## 20191128 try to pca plot for 8 province with Methylation

##20200323try EWAS y1-1　
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS")
load("dataGeo.RData")
load("dataMethy.RData")
gc()
blogit <- function(x)
{
# logistic model
fit <- glm(dataGeo$y1 ~ x + dataGeo$age, family=binomial())
s <- summary(fit)
blogitRe <- coef(s)[2, c(1,2,3,4)]
return(blogitRe)
}
##error from x length differ from dataGeo
 which(colnames(usedMethy) %in% dataGeo$IID)
 dataMethy <- usedMethy[,which(colnames(usedMethy) %in% dataGeo$IID)]
 ## keep dataGeo and dataMethy dim equle

 dataGeo$y1 <- dataGeo$y1-1
y1_dataMethy_bligitResult <- apply(dataMethy, 1, blogit)
write.table(y1_dataMethy_bligitResult, file = "y1_dataMethy_bligit.Result", quote = F, row.names = F)

library(qqman)
png("y1_dataMethy_bligit_age.Result.qqplot.png")
qq(y1_dataMethy_bligitResult[4, ], main = "Q-Q plot of y1_dataMethy")
dev.off()
load("/thinker/storage/org/liufanGroup/fanxiu/geography/190809Geography/usedAno.RData")
gc()
methyID_logitResult <- data.frame(colnames(y1_dataMethy_bligitResult),t(y1_dataMethy_bligitResult))
usedManhattan <- merge(usedAno, methyID_logitResult, by.x ="probID", by.y = "colnames.y1_dataMethy_bligitResult.")
usedManhattanplot <- data.frame(usedManhattan$probID, as.integer(usedManhattan$CHR), usedManhattan$BP, usedManhattan$Pr...z..)
names(usedManhattanplot) <- c("cpg", "CHR", "BP", "P")
gc()
png("y1_dataMethy_bligit_age.Result.manhattan.png")
manhattan(usedManhattanplot, col = c("blue4", "red3"), cex = 0.6, cex.axis = 0.9, ymax = 12, main = "manhattanplot．y1_dataMethy")
dev.off()
##20200323try EWAS y1-2 conditional cg16529197
which.min(y1_dataMethy_bligitResult[4,])
con1 <- dataMethy[which.min(y1_dataMethy_bligitResult[4,]),]
##
gc()
blogit <- function(x)
{
# logistic model
fit <- glm(dataGeo$y1 ~ x + dataGeo$age + dataGeo$con1, family=binomial())
s <- summary(fit)
blogitRe <- coef(s)[2, c(1,2,3,4)]
return(blogitRe)
}
##error from x length differ from dataGeo
 #which(colnames(usedMethy) %in% dataGeo$IID)
# dataMethy <- usedMethy[,which(colnames(usedMethy) %in% dataGeo$IID)]
 ## keep dataGeo and dataMethy dim equle

y1_dataMethy_bligitResult_2 <- apply(dataMethy, 1, blogit)
write.table(y1_dataMethy_bligitResult_2, file = "y1_dataMethy_bligit.Result_2", quote = F, row.names = F)
library(data.table)
y1_dataMethy_bligitResult_2 <- fread("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS/y1_dataMethy_bligit.Result_2", header = TRUE)
which.min(y1_dataMethy_bligitResult_2[4,])
##!!!!find min(result2)= min(result1) cg16329197
## try 2 y1~x+age+conCPG1
library(data.table)
load("dataGeo.RData")
load("dataMethy.RData")
gc()
y1_dataMethy_bligitResult_1 <- fread("y1_dataMethy_bligit_age.Result_1", header = TRUE)
conCPG1 <- dataMethy[which.min(y1_dataMethy_bligitResult_1[4,]),]
dataGeo$conCPG1 <- conCPG1
blogit <- function(x)
{
# logistic model
fit <- glm(dataGeo$y1 ~ x + dataGeo$age + dataGeo$conCPG1, family=binomial())
s <- summary(fit)
blogitRe <- coef(s)[2, c(1,2,3,4)]
return(blogitRe)
}
y1_dataMethy_bligit_age.Result_2 <- apply(dataMethy, 1, blogit)
write.table(y1_dataMethy_bligit_age.Result_2, file = "y1_dataMethy_bligit_age.Result_2", quote = F, row.names = F)
which.min(y1_dataMethy_bligit_age.Result_2[4,])
##!!min(result2)=cg16329197
fit1 <- glm(dataGeo$y1~dataMethy[588749,]+dataGeo$age+dataGeo$conCPG1, family = binomial())
testMethy <- data.frame(t(dataMethy[588749:588755, ]))
lapply(testMethy, )
formular <- "dataGeo$y1 ~ x + dataGeo$age + dataGeo$conCPG1"
lapply(testMethy,function(x){
                 coef(summary(glm(dataGeo$y1 ~ x + dataGeo$age + dataGeo$conCPG1, family=binomial())))[2, c(1,2,3,4)]
		})
#!x               40.629610   3.024908  13.432  < 2e-16 ***
#dataGeo$age      0.023984   0.007758   3.091  0.00199 **
#dataGeo$conCPG1        NA         NA      NA       NA
		
lapply(testMethy[,1:2],function(x){
                 summary(glm(dataGeo$y1 ~ dataGeo$age + dataGeo$conCPG1 + x, family=binomial()))
		})
#dataGeo$age      0.023984   0.007758   3.091  0.00199 **
#dataGeo$conCPG1 40.629610   3.024908  13.432  < 2e-16 ***
#@x                      NA         NA      NA       NA
y1_dataMethy_bligit_age.Result_2[,588749] <- NA
y1_dataMethy_bligit_age.Result_2[,482197]
y1_dataMethy_bligit_age.Result_2 <- lapply(dataMethy,function(x){
                 coef(summary(glm(dataGeo$y1 ~ dataGeo$age + dataGeo$conCPG1 + x, family=binomial())))[2, c(1,2,3,4)]
		})
write.table(y1_dataMethy_bligit_age.Result_2, file = "y1_dataMethy_bligit_age.Result_2", quote = F, row.names = F)
##!lapply data.frame
dataMethy <- data.frame(t(dataMethy))
testMethy <- data.frame(t(dataMethy[588749:588755, ]))
formu <- formula(paste(y1,"~ x"))
lapply(testMethy,function(x){
                formu <- formula(paste(y1,"~ x"))
                 coef(summary(glm(formu, family=binomial())))[2, c(1,2,3,4)]
		})
		         formu <- formula(paste(dataGeo$y1,"~ x"))
				 tmp <- cbind(x, dataGeo)
logitfun <- function(x)
        {
            tmp <- cbind(x, dataGeo)
            summary(glm(fo, data = tmp, family = binomial(link = "logit"), control = list(maxit = logitmaxit)))$coef[2,]
        }
y <- factor(dataGeo$y1)
	fo <- formula(paste(y,"~ x"))		
 logitfun <- function(x)
        {
            tmp <- cbind(dataGeo,x)
            summary(glm(fo, data = tmp, family = binomial(link = "logit"), control = list(maxit = logitmaxit)))$coef[2,]
        }
        ewasresult <- data.frame(t(lapply(testMethy, logitfun)))
		
## plot result  y1_ 2 
library(qqman)
png("y1_dataMethy_bligit_age.Result_2.qqplot.png")
qq(as.numeric(y1_dataMethy_bligit_age.Result_2[4, ]), main = "Q-Q plot of y1_dataMethy_2")
dev.off()
load("/thinker/storage/org/liufanGroup/fanxiu/geography/190809Geography/usedAno.RData")
gc()
methyID_logitResult <- data.frame(colnames(y1_dataMethy_bligit_age.Result_2),t(y1_dataMethy_bligit_age.Result_2))
usedManhattan <- merge(usedAno, methyID_logitResult, by.x ="probID", by.y = "colnames.y1_dataMethy_bligit_age.Result_2.")
usedManhattanplot <- data.frame(usedManhattan$probID, as.integer(usedManhattan$CHR), usedManhattan$BP, usedManhattan$Pr...z..)
names(usedManhattanplot) <- c("cpg", "CHR", "BP", "P")
gc()
png("y1_dataMethy_bligit_age.Result_2.manhattan.png")
manhattan(usedManhattanplot, col = c("blue4", "red3"), cex = 0.6, cex.axis = 0.9, ymax = 12, main = "manhattanplot．y1_dataMethy_2")
dev.off()
## 20200324 rearrange data geo male 990
data990 <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS/GABmale990Methy.txt", header = TRUE)
dataProvin <- na.omit(data990)
dataProvin$y1 <- dataProvin$y1-1
dataProvin$y1 <- factor(dataProvin$y1-1)
dataProvin$y2 <- factor(dataProvin$y2-1)
dataProvin$y3 <- factor(dataProvin$y3-1)
dataProvin$y4 <- factor(dataProvin$y4-1)
dataProvin$y5 <- factor(dataProvin$y5-1)
dataProvin$y6 <- factor(dataProvin$y6-1)
dataProvin$y7 <- factor(dataProvin$y7-1)
dataProvin$y8 <- factor(dataProvin$y8-1)
save(dataProvin, file = "dataProvin.RData") # 972  14

load("/thinker/storage/org/liufanGroup/fanxiu/geography/190809Geography/usedMethy.RData")
dataCpG <- usedMethy[,which(colnames(usedMethy) %in% dataProvin$IID)] #723730    972
save(dataCpG, file = "dataCpG.RData")	

##20200329
#nohup /LiuFanBin/R-3.6.1 < /thinker/storage/org/liufanGroup/fanxiu/ancestry/GAB_1000G/2509GAB1000G/ConditionalGWAS.Merge/LRLogis111SNP/leave1out.LR.R --no-save > /thinker/storage/org/liufanGroup/fanxiu/ancestry/GAB_1000G/2509GAB1000G/ConditionalGWAS.Merge/LRLogis111SNP/leave1out.LR.R.log &
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1')) 
library(qqman)
######################################################################################
##20200409
##try predict province by languaugeFamily 111SNP
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince
R
prov8 <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/province8.GAB2509.1000G400.snp111.txt", header = TRUE)
prov8$y <- factor(prov8$y)
save(prov8, file = "prov8.RData")
usedPro <- prov8[,-c(1:4)]
dim(na.omit(usedPro)) #2838*112
table(usedPro$y)
usedPro <- na.omit(usedPro)
library(dplyr)
library(caret)
library("nnet")
library(pROC)
logis1 <- multinom(y~., data = usedPro, MaxNWts = 2000)
pre <- predict(logis1)
per_logis1 <-  confusionMatrix(factor(pre),usedPro$y)
per_logis1
prov8$y1 <- ifelse(prov8$y==1,2,1) 
prov8$y2 <- ifelse(prov8$y==2,2,1) 
prov8$y3 <- ifelse(prov8$y==3,2,1) 
prov8$y4 <- ifelse(prov8$y==4,2,1) 
prov8$y5 <- ifelse(prov8$y==5,2,1) 
prov8$y6 <- ifelse(prov8$y==6,2,1) 
prov8$y7 <- ifelse(prov8$y==7,2,1) 
prov8$y8 <- ifelse(prov8$y==8,2,1) 
prov8$y9 <- ifelse(prov8$y==9,2,1) 
prov8$y10 <- ifelse(prov8$y==10,2,1) 
prov8$y11 <- ifelse(prov8$y==11,2,1) 
prov8$y12 <- ifelse(prov8$y==12,2,1)
library(dplyr) 
prov8_2 <- select(prov8, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, everything())
write.table(prov8_2, file = "prov8.y1.y8.snp11", quote = FALSE, row.names = FALSE)
usedPro <- prov8_2[,-c(13:16)]
usedPro <- na.omit(usedPro)
dim(usedPro)# 2838*124
prob = predict(logis1, type='probs') #2838*12
auc <- c(NA)
for (k in 1:12){
          #calculate 12 auc for each province 
		  roc_result <- roc(usedPro[, k]-1, prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
		  auc
##20200410 leaave one out precict by snp111
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince
R
load("prov8.RData")
prov8$y1 <- ifelse(prov8$y==1,2,1) 
prov8$y2 <- ifelse(prov8$y==2,2,1) 
prov8$y3 <- ifelse(prov8$y==3,2,1) 
prov8$y4 <- ifelse(prov8$y==4,2,1) 
prov8$y5 <- ifelse(prov8$y==5,2,1) 
prov8$y6 <- ifelse(prov8$y==6,2,1) 
prov8$y7 <- ifelse(prov8$y==7,2,1) 
prov8$y8 <- ifelse(prov8$y==8,2,1) 
prov8$y9 <- ifelse(prov8$y==9,2,1) 
prov8$y10 <- ifelse(prov8$y==10,2,1) 
prov8$y11 <- ifelse(prov8$y==11,2,1) 
prov8$y12 <- ifelse(prov8$y==12,2,1)
library(dplyr) 
prov8 <- select(prov8, IID, sex, population, place, y, y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, everything())
save(prov8, file = "prov8.RData")
##leave1out
##20200410 multinom logistic regression ( !!!!!!!!!!!!!Add  MaxNWts = 2000 )
## GAB2509+1000G400 prov8+4  snp111
#nohup /LiuFanBin/R-3.6.1 < /thinker/storage/org/liufanGroup/fanxiu/ancestry/GAB_1000G/2509GAB1000G/ConditionalGWAS.Merge/leave1outAUC/auc.maxit1000.GAB_1000G.R --no-save > /thinker/storage/org/liufanGroup/fanxiu/ancestry/GAB_1000G/2509GAB1000G/ConditionalGWAS.Merge/leave1outAUC/auc.maxit1000.GAB_1000G.R.log &
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1')) 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince")
library(dplyr)
library("nnet")
library(pROC)
load("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/prov8.RData")
fitdata <- na.omit(prov8[,-c(1:4)])
snpname <- colnames(fitdata)[14:ncol(fitdata)]
yname <- as.character(paste("y", c(1:12), sep = ""))
auc <- matrix(c(NA), length(yname))
pre <- c(NA)
     fo <- as.formula(paste("y ~ ", paste(snpname[1:length(snpname)], collapse = "+")))
     prob <- matrix(NA, nrow(fitdata), length(yname))
     for(j in 1:nrow(fitdata)){
          #leave one for validation others for training
          validation<-fitdata[j,]
          training<- fitdata[-j,]
		  # fit model 
          logis1 <- multinom(fo, data = training, trace = FALSE, MaxNWts = 2000)
		  # predict 12 probility 
		  prob[j,] = predict(logis1, newdata = validation, type='probs')
		  pre[j] <- predict(logis1, newdata = validation)
          write.table(prob, file= "MaxNWts2000.pre", row.names = FALSE, quote = FALSE)
		  write.table(pre, file= "MaxNWts2000.pro", row.names = FALSE, quote = FALSE)
		  }
     for (k in 1:length(yname)){
          #calculate 12 auc for each language family class 
		  roc_result <- roc(fitdata[, k+1], prob[,k], auc=T, quiet = TRUE)
          auc[k] <- roc_result$auc
          } 
write.table(auc, file = "AUC.MaxNWts2000.leave1out.prov8.SNP", row.name = F, quote = F)
pro <- read.table("MaxNWts2000.pro", header = TRUE)
pre <- read.table("MaxNWts2000.pre", header = TRUE)
per_logis1 <-  confusionMatrix(factor(pro$x),factor(fitdata$y))
per_logis1
pre <- as.matrix(pre)
auc <- c(NA)
for (k in 1:12){
      roc_result <- roc(fitdata[,k+1], pre[,k], auc = TRUE, quiet = TRUE)
      auc[k] <- roc_result$auc
}
auc
## MeQTL predict province
cd /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result
head /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/cisSigRlt
head /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/select111SNP.list
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince
awk 'FNR==NR{a[$1]=$0;next}{ if(a[$1]) {print $0, a[$1]} }' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/select111SNP.list /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/cisSigRlt > cismatch.snp111.cpg
awk 'FNR==NR{a[$1]=$0;next}{ if(a[$1]) {print $0, a[$1]} }' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/select111SNP.list /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/transSigRlt > transmatch.snp111.cpg
## after select Cpg 360
## predict8 province 
R 
load("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS/dataCpG.RData")
load("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/conditionalEWAS/dataProvin.RData")
dim(dataCpG)  ## 723730    972
cpg360 <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/CpG360.list", header = TRUE)
dim(cpg360)  ## 360   1
cpg360 <- cpg360$DNAm
selectCpG <- dataCpG[rownames(dataCpG) %in% cpg360,]
dim(selectCpG)  ## 360 972
IID <- colnames(selectCpG)
selectCpG <- data.frame(IID,t(selectCpG))
dim(selectCpG)  ##972 361
prov8.CpG360 <- merge(dataProvin, selectCpG, by = "IID")
dim(prov8.CpG360)  ##972 374
write.table(prov8.CpG360, file = "prov8.CpG360.GABNna972", quote = FALSE, row.names = FALSE)
save(prov8.CpG360, file = "prov8.CpG360.RData")
## try multilogis
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis")
load("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis/prov8.CpG360.RData")
data <- prov8.CpG360[,-c(1:5)]
cpgname <- colnames(data)[10:ncol(data)]
     fo <- as.formula(paste("province ~ ", paste(cpgname[1:length(cpgname)], collapse = "+")))
library(nnet)
logis1 <- multinom(fo, data = prov8.CpG360, trace = FALSE,MaxNWts = 3000)
pre <- predict(logis1)
per_logis1 <-  confusionMatrix(factor(pre),factor(data$province))
prob = predict(logis1, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(data[, k+1], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
## cpg360 predict prov 8 all correct GAB972
## test leave one out 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis")
library(dplyr)
library("nnet")
library(pROC)
load("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis/prov8.CpG360.RData")
fitdata <- na.omit(prov8.CpG360[,-c(1:5)])
cpgname <- colnames(fitdata)[10:ncol(fitdata)]
yname <- as.character(paste("y", c(1:8), sep = ""))
auc <- c(NA)
pre <- c(NA)
     fo <- as.formula(paste("province ~ ", paste(cpgname[1:length(cpgname)], collapse = "+")))
     prob <- matrix(NA, nrow(fitdata), length(yname))
     for(j in 1:nrow(fitdata)){
          #leave one for validation others for training
          validation<-fitdata[j,]
          training<- fitdata[-j,]
		  # fit model 
          logis1 <- multinom(fo, data = training, trace = FALSE, MaxNWts = 3000)
		  # predict 8 probility 
		  prob[j,] = predict(logis1, newdata = validation, type='probs')
		  pre[j] <- predict(logis1, newdata = validation)
          write.table(prob, file= "MaxNWts3000.cpg360.pre", row.names = FALSE, quote = FALSE)
		  write.table(pre, file= "MaxNWts3000.cpg360.pro", row.names = FALSE, quote = FALSE)
		  }
     for (k in 1:length(yname)){
          #calculate 8 auc for each language family class 
		  roc_result <- roc(fitdata[, k+1], prob[,k], auc=T, quiet = TRUE)
          auc[k] <- roc_result$auc
          } 
write.table(auc, file = "AUC.MaxNWts3000.leave1out.prov8.cpg360", row.name = F, quote = F)

##20200411 ## manage leave one out CpG360 performance
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis")
pro <- read.table("MaxNWts3000.cpg360.pro", header = TRUE)
pre <- read.table("MaxNWts3000.cpg360.pre", header = TRUE)
per_logis1 <-  confusionMatrix(factor(pro$x),factor(fitdata$province))
per_logis1
pre <- as.matrix(pre)
#################################################
##20200411
## try mean(auc) 8 -2 select CpG360
# cpg360  in
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis")
load("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/pro8.cpg360.multilogis/prov8.CpG360.RData") 
fitdata <- na.omit(prov8.CpG360[,-c(1:5)])
cpgname <- colnames(fitdata)[10:ncol(fitdata)]
meanAUC <- c()
meanAUC[1] <- 0.5
yname <- "y1"
##这里应该添加对y1二分类多重回归的P值排序，按照由小到大进入下面的模型，
fo <- as.formula(paste(yname," ~ ", paste(cpgname[1:length(cpgname)], collapse = "+")))
fit1 <- glm(fo, data = fitdata, family=binomial())
s <- summary(fit1)
blogitRe <- coef(s)[-1, c(1,2,3,4)]
write.table(blogitRe, file = "result.y1.blogit.cpg360", quote = FALSE, row.names = FALSE)
##!!!
#Warning messages:
#1: glm.fit: algorithm did not converge
#2: glm.fit: fitted probabilities numerically 0 or 1 occurred
#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
#(Intercept)  8.493e+01  9.314e+06       0        1
#cg14003978  -1.352e+01  7.406e+05       0        1
##所以不排序直接挑选对于每一个省份
library(dplyr)
library(pROC)
for (k in 1:length(cpgname)){
     fo <- as.formula(paste(yname," ~ ", paste(cpgname[1:k], collapse = "+")))
     auc <- c()
     for (j in 1:100){
          #100 times 8-2 cross validation -> mean(AUC)
          N <- nrow(fitdata)
          folds <- 5
          holdout <- split(sample(1:N), 1:folds)
          holdout %>% unlist() %>% length() == N
          fit <- glm(fo, data = fitdata[-holdout$`1`,], family=binomial())
          probability <- predict(fit, newdata = fitdata[holdout$`1`,])
          roc_result <- roc(factor(fitdata[holdout$`1`,]$y1), probability, auc=T, quiet = TRUE)
          auc[j] <- roc_result$auc
          }
     meanAUC[k+1] <- mean(auc)
     if(meanAUC[k+1] < meanAUC[k]){ 
	    break
       }
	   k
}
selectCpG <- cpgname[1:k]
write.table(selectCpG, file = "selectCpG.meanAUC.province.y1", quote = FALSE, rownames = FALSE)

meanAUC <- meanAUC[-1]
png("/thinker/storage/org/liufanGroup/fanxiu/ancestry/try5_intersect3qc/multilogit/ancestry.y3Minor.10snp.AUC.png")
plot(meanAUC, main="province.y1.cpg360.AUC", xlab="number.cpg ", ylab="mean.AUC ", pch=19, type = "l", xlim=c(0,10), ylim=c(0,1))
x= 1:length(meanAUC)
for (i in 1:length(meanAUC)){
     text(x[i],meanAUC[i],as.character(round(meanAUC[i],3)), cex=1.2)
}
dev.off()
y3selectSNP <- snp1000y3[pmatch(snpin[1:k],snp1000y3$SNP),]
write.table(y3select10SNP, file= "/thinker/storage/org/liufanGroup/fanxiu/ancestry/try5_intersect3qc/multilogit/y3.minor.select.SNP", quote = FALSE, row.names = FALSE)
###############################
##20200413 try Cpg in according AUC增量由大到小
# try 1SNP ~ 1CpG而不是用多个CpG
datacpg <- fitdata[, 10:ncol(fitdata)]
corcpg <- cor(datacpg)^2
color <- gray(100:0/100)
tiff("cor2.cpg360.tiff")
pheatmap(corcpg, cluster_rows= F, cluster_cols = F, col = color, labels_row = "", labels_col="")
dev.off()
color <- heat.colors(100:0/100)
pheatmap(corcpg[1：10， 1：10], cluster_rows= F, cluster_cols = F, col = color)
##
a <- which.max(corcpg)
a <- apply(corcpg,2,which.max)	

######################################################################
##20200419## JXX LXX conditional EWAS try 2
## data prepare 
/LiuFanBin/R-3.6.1
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1')) 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang")
data <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/GABMethydata.csv", header = T)
dim(data) # 990*13
id <- gsub("-", "_", data$Me.id)
data$Me.id <- id
write.csv(data, file = "GABmale985Methy.csv", quote = F, row.names = F)
datapro <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/GABmale985Methy.csv", header = TRUE)
datacov <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/defcov.csv", header = TRUE)
dataphe <- merge(datapro, datacov, by = "Me.id")
dataphe$y1 <- factor(dataphe$y1-1)
dataphe$y2 <- factor(dataphe$y2-1)
dataphe$y3 <- factor(dataphe$y3-1)
dataphe$y4 <- factor(dataphe$y4-1)
dataphe$y5 <- factor(dataphe$y5-1)
dataphe$y6 <- factor(dataphe$y6-1)
dataphe$y7 <- factor(dataphe$y7-1)
dataphe$y8 <- factor(dataphe$y8-1)
write.csv(dataphe, file = "usedCon2.GABmale985Methy.csv", quote = F, row.names = F)
#ewasFun 
	  ewasFun <- function(usenormbetadata = FALSE, phedata, selphe, covariates = c("B","NK","CD4T","CD8T","Mono","Neutro","Eosino","PC1","PC2","PC3","PC4","PC5"), type = "linear", isbetadependent = FALSE, savepath = ".", logitmaxit = 25, chunknum = 1, plotmq = TRUE)
{
    cat("\n######")
    cat("\n")
    if(length(selphe) != 1){
        stop("selphe must be a single string")
    }
    if(!dir.exists(savepath)){
        message("\n###\nMake directory: ", savepath)
        dir.create(file.path(savepath), recursive = TRUE)
    }
    message("\n###\nLoad the packages and data\n") #load data
    for(pack in c("parallel", "data.table", "dplyr", "qqman", "plyr")){
        suppressPackageStartupMessages(library(pack, character.only = TRUE))
    }
    phedata <- read.csv(file = phedata, as.is = TRUE)
    defcov <- read.csv("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/defcov.csv", as.is = TRUE) #default covariates
    samephe <- names(defcov) %in% names(phedata)
    samephe[1] <- FALSE
    if(sum(samephe) > 0){
        message("Attention: your phedata has the following same variates with the default loading covariates file, and the ones in your phedata will be used instead of the default:\n", paste(names(defcov)[samephe], collapse = ", "), "\n")
        defcov <- subset(defcov, select = !samephe)
    }
    phedata <- merge(defcov, phedata, by = "Me.id") #merge user supplied phedata
    if(!all(c(selphe, covariates) %in% names(phedata)))
    {
        stop("selphe and covariates must be contained in the defcov.csv or phedata")
    }
    load("/liufanGroup/liuxx/DNAm_data/genoinfo.RData") #genomic information
    if(usenormbetadata){
        load("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/norm_order_betadata.RData") #normalized and ordered beta values data
        betadata <- norm_order_betadata
        rm(norm_order_betadata)
        gcresult <- gc()
    } else{
        load("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/betadata.RData") #primiry beta values data
    }
    message("###\nPreprocess data\n") #select and match data
    phedata <- subset(phedata, select = c("Me.id", selphe, covariates))
    isuncompleterow <- !complete.cases(phedata)
    if(any(isuncompleterow)){
        if(sum(isuncompleterow) < 6){
            message("The individuals: ", paste(phedata$Me.id[isuncompleterow], collapse = ", "), " have missing values in the phenotype you selected or covariates, so they will be discarded")
        } else{
            message(message("The individuals: ", paste(phedata$Me.id[isuncompleterow][1:5], collapse = ", "), "... (total number: ", sum(isuncompleterow), ") have missing value in the phenotype you selected or covariates, so they will be discarded"))
        }
        writeLines(phedata$Me.id[isuncompleterow], paste0(savepath, "/", selphe, "orCovariates_IndividualsHaveMissingValue.txt"))
        phedata <- na.omit(phedata)
    }
    nothavebetadata <- phedata$Me.id[!(phedata$Me.id %in% rownames(betadata))]
    if(length(nothavebetadata) > 0){
        if(length(nothavebetadata) < 6){
            message("The individuals: ", paste(nothavebetadata, collapse = ", "), " have not data in the default loading methylation beta data, so they will be discarded")
        } else{
            message(message("The individuals: ", paste(nothavebetadata, collapse = ", "), "... (total number: ", length(nothavebetadata), ") have not data in the default loading methylation beta data, so they will be discarded"))
        }
        writeLines(nothavebetadata, paste0(savepath, "/", "IndividualsHaveNotBetaData.txt"))
    }
    betadata <- betadata[match(phedata$Me.id, rownames(betadata), nomatch = 0), ] #select and arrange betadata
    phedata <- join(data.frame(Me.id = rownames(betadata)), phedata, by = "Me.id") #select and arrange phedata in the order of id in betadata
    phedata <- subset(phedata, select = c(selphe, covariates)) #take selphe to the first column and remove Me.id.
    message("###\nStart analysis\n") #start regression
    if(isbetadependent){
        message("Use each beta value as the dependent variable\n")
        fo <- formula("x ~ .") #specify each beta value as the dependent variable and the phenotype as the predictor of interest
    } else{
        message("Use ", selphe, " as the dependent variable\n")
        fo <- formula(paste(selphe, "~ .")) #use phenotype as the dependent variable
    }
    linearfun <- function(x) #linear regression
    {
        tmp <- cbind(x, phedata)
        smodel <- summary(lm(fo, data = tmp))
        c(smodel$coef[2, ], R.squared = smodel$r.squared)
    }
    logitfun <- function(x) #logit regression
    {
        tmp <- cbind(x, phedata)
        summary(glm(fo, data = tmp, family = binomial(link = "logit"), control = list(maxit = logitmaxit)))$coef[2, ]
    }
    slicecol <- ceiling(seq.int(1, ncol(betadata) + 1, length.out = chunknum + 1)) #pre-split the data into chunknum chunks
    my.env <- new.env() #save the split data in an environment
    for(i in seq_len(chunknum)){
        assign(paste0("data", i), betadata[, slicecol[i]:(slicecol[i+1] - 1)], envir = my.env)
    }
    message("Use ", type, " regression\n")
    applyMyFun <- function(idx, env){ #function that takes indexes, then extracts the data from the environment
        eval(parse(text = paste0("tmpdata <- env$", ls(env)[idx])))
        data.frame(t(sapply(tmpdata, get(paste0(type, "fun")))))
    }
    index <- 1:chunknum
    names(index) <- 1:chunknum
    ewasresult <- mclapply(index, applyMyFun, env = my.env, mc.cores = chunknum)
    ewasresult <- do.call("rbind", ewasresult)
    message("###\nOrganize result\n") #rename and arrange the ewasresult
    names(ewasresult)[1:4] <- c("Beta", "SE", "STAT", "P")
    ewasresult$probid <- gsub(".*\\.", "", rownames(ewasresult))
    rownames(ewasresult) <- NULL
    ewasresult$FDR.P <- p.adjust(ewasresult$P, method = "fdr")
    rownames(genoinfo) <- genoinfo[, 1]
    ewasresult <- data.frame(ewasresult, genoinfo[ewasresult$probid, 3:4])
    rownames(ewasresult) <- NULL
    ewasresult <- select(ewasresult, probid, everything())
    ewasresult$N <- nrow(phedata)
    ewasresult <- arrange(ewasresult, P)
    message("Add genomic control\n") #genomic control--Correction for Population Stratification
    ewasresult$gc.P <- NA
    n = nrow(betadata)
    tSTAT <- ewasresult$STAT[!is.na(ewasresult$P)]
    tN <- ewasresult$N[!is.na(ewasresult$P)]
    expected_stat <- rt(length(tSTAT), n)
    genomic_control <- median(abs(tSTAT)) / median(abs(expected_stat))
    tSTAT <- tSTAT / genomic_control
    ewasresult$gc.P[seq_along(tSTAT)][tSTAT < 0] <- pt(tSTAT[tSTAT < 0], tN[tSTAT < 0] - 1) * 2
    ewasresult$gc.P[seq_along(tSTAT)][tSTAT > 0] <-(1 - pt(tSTAT[tSTAT > 0], tN[tSTAT > 0] - 1)) * 2
    message("###\nSave EWAS Results in", savepath, "/ewasresult_", selphe, ".csv...\n") #EWAS results
    fwrite(ewasresult, file = paste0(savepath, "/ewasresult_", selphe, ".csv"), quote = FALSE)
    message("Done\n")
    if(plotmq){
        message("###\nPlot manhattan and qq plots without genomic control...\n") #manhattan and qq plots without genomic control
        index_Piszero <- which(ewasresult$P == 0)
        if(length(index_Piszero) != 0){
            if(length(index_Piszero) <= 5){
                message("Attention: ", paste(ewasresult$probid[index_Piszero], collapse = ", "), "  P = 0\n")
            } else{
                message("Attention: ", paste(ewasresult$probid[index_Piszero][1:5], collapse = ", "), "... (total number: ",length(index_Piszero), ")  P = 0\n")
            }
            message("To draw the manhattan plot, change the P = 0 to P = 9.881313e-324 temporarily, other results are not affected\n")
            ewasresult$P[index_Piszero] <- 9.881313e-324
            writeLines(ewasresult$probid[index_Piszero], paste0(savepath, "/", selphe, "_PValueisZero.txt"))
        }
        tiff(filename = paste0(savepath, "/manhattan_plot_", selphe, ".tif"), bg = "white", width = 3000, height = 1500, res = 240)
        manhattan(x = ewasresult, snp = "probid", cex = 0.6, cex.axis = 0.6, col = c("red", "blue"), suggestiveline = -log10(1e-06), main = paste0("manhattan_plot_", selphe))
        dev.off()
        ewasresult$P[index_Piszero] <- 0
        tiff(filename = paste0(savepath, "/qq_plot_", selphe, ".tif"), bg = "white", width = 1500, height = 1500, res = 240)
        qq(ewasresult$P, main = paste0("qq_plot_", selphe))
        dev.off()
        message("Done\n")
        message("###\nPlot manhattan and qq plots with genomic control...\n") #manhattan and qq plots after genomic control
        index_gc.Piszero <- which(ewasresult$gc.P == 0)
        if(length(index_gc.Piszero) != 0){
            if(length(index_gc.Piszero) <= 5){
                message("Attention: ", paste(ewasresult$probid[index_gc.Piszero], collapse = ", "), "  gc.P = 0\n")
            } else{
                message("Attention: ", paste(ewasresult$probid[index_gc.Piszero][1:5], collapse = ", "), "... (total number: ",length(index_gc.Piszero), ")  gc.P = 0\n")
            }
            message("To draw the manhattan plot, change the gc.P = 0 to gc.P = 9.881313e-324 temporarily, other results are not affected\n")
            ewasresult$gc.P[index_gc.Piszero] <- 9.881313e-324
            writeLines(ewasresult$probid[index_gc.Piszero], paste0(savepath, "/", selphe, "_gc.PValueisZero.txt"))
        }
        tiff(filename = paste0(savepath, "/manhattan_gc_plot_", selphe, ".tif"), bg = "white", width = 3000, height = 1500, res = 240)
        manhattan(x = ewasresult, p = "gc.P", snp = "probid", cex = 0.6, cex.axis = 0.6, col = c("red", "blue"), suggestiveline = -log10(1e-06), main = paste0("manhattan_gc_plot_", selphe))
        dev.off()
        ewasresult$gc.P[index_gc.Piszero] <- 0
        tiff(filename = paste0(savepath, "/qq_gc_plot_", selphe, ".tif"), bg = "white", width = 1500, height = 1500, res = 240)
        qq(ewasresult$gc.P, main = paste0("qq_gc_plot_", selphe))
        dev.off()
        message("Done\n")
    }
    message("###\n", selphe, " EWAS (", type, ") results are saved in ", savepath, "\n\n######\n")
}
##test y1_1
ewasFun(usenormbetadata = TRUE, phedata = "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/usedCon2.GABmale985Methy.csv", selphe = "y1", covariates = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age"), type = "logit", isbetadependent = FALSE, savepath = "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/ewasFun_output", logitmaxit = 25, chunknum = 1, plotmq = TRUE)
############################################################################################
##try to select
/LiuFanBin/R-3.6.1
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1'))  
library(data.table())
workpath <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/"
setwd(workpath)
# selectCpG <- NA
yname = "y1"
covnames = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age")
resultname <- "province.GABmale.conditional.ewas1e-6.y1"
selectCPG <- matrix(NA, 20,8)
for(i in 1: 20){
    if(i ==1){
	    phedata <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/usedCon2.GABmale985Methy.csv"
        savepath <- paste0(workpath, "y1CellAgeresult_",i)
        ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
        #filter P < 1e-6 save
	    tmp <- paste0("awk -F, '$5<1e-6 {print $1,$7,$8,$5}' ",savepath, "/ewasresult_",yname,".csv"," > ",workpath,resultname,"_",i)
		system(tmp)
		# plotManhattan resultfile= resultname_i
	    resultfile <- paste0(resultname,"_",i)
	    dataManhattan <- fread(resultfile, header = FALSE)
		# if no CPG P <1e-6 break
	    if(nrow(dataManhattan)==0){
		           message("no CPG P <1e-6\n")
	               break
		               }
	    names(dataManhattan) <- c("CPG", "CHR", "BP", "P")
	    # select top snpfile = resultname_i.selectSNP_i(snplist).
	   selectCPG[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(dataManhattan$P)),c(1:4)])
        #save
	   CPGfile <- paste0(resultfile,".selectCPG_",i)
	   write.table(selectCPG, file = CPGfile, quote = FALSE, row.names = FALSE)
	   condname <- as.character(selectCPG[i,1])
       rm(dataManhattan)
       gc()
	   ##get conditional CPG value
	   load("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/norm_order_betadata.RData")
       condvalue <- data.frame(rownames(norm_order_betadata),norm_order_betadata[,which(colnames(norm_order_betadata)==condname)])
       names(condvalue) <- c("Me.id", condname)
	   phevalue <- read.csv(phedata, header = TRUE)
	   phecondata <- merge(phevalue, condvalue, by = "Me.id")
	   # save next pheconddata
	   pheconfile <- paste0(resultfile,".phecondata_",i)
       write.csv(phecondata, file = pheconfile, quote = FALSE, row.names=FALSE)
	}else{
	    phedata <- paste0(workpath,pheconfile)
        savepath <- paste0(workpath, "y1CellAgeresult_",i)
		covnames <- c(covnames,condname)
        ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
        #filter P < 1e-6 save
	    tmp <- paste0("awk -F, '$5<1e-6 {print $1,$7,$8,$5}' ",savepath, "/ewasresult_",yname,".csv"," > ",workpath,resultname,"_",i)
		system(tmp)
		# plotManhattan resultfile= resultname_i
	    resultfile <- paste0(resultname,"_",i)
	    dataManhattan <- fread(resultfile, header = FALSE)
		# if no CPG P <1e-6 break
	    if(nrow(dataManhattan)==0){
		           message("no CPG P <1e-6\n")
	               break
		               }
	    names(dataManhattan) <- c("CPG", "CHR", "BP", "P")
	    # select top snpfile = resultname_i.selectSNP_i(snplist).
	    selectCPG[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(dataManhattan$P)),c(1:4)])
        #save
	   CPGfile <- paste0(resultfile,".selectCPG_",i)
	   write.table(selectCPG, file = CPGfile, quote = FALSE, row.names = FALSE)
	   condname <- as.character(selectCPG[i,1])
       rm(dataManhattan)
       gc()
	   ##get conditional CPG value
	   load("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/norm_order_betadata.RData")
       condvalue <- data.frame(rownames(norm_order_betadata),norm_order_betadata[,which(colnames(norm_order_betadata)==condname)])
       names(condvalue) <- c("Me.id", condname)
	   phevalue <- read.csv(phedata, header = TRUE)
	   phecondata <- merge(phevalue, condvalue, by = "Me.id")
	   # save next pheconddata
	   pheconfile <- paste0(resultfile,".phecondata_",i)
       write.table(phecondata, file = pheconfile, quote = FALSE, row.names=FALSE)
	}
	}
	############################################################
	##20200422 gm132-1 restart
	#test neimeng y2
	> warnings()
Warning messages:
1: glm.fit: fitted probabilities numerically 0 or 1 occurred
2: glm.fit: fitted probabilities numerically 0 or 1 occurred
gc()
 used (Mb) gc trigger    (Mb)   max used    (Mb)
Ncells  612371 32.8    6505373   347.5    8131716   434.3
Vcells 3254505 24.9 1402152811 10697.6 1752665689 13371.8
#####!!!!!!!!!!!!!!!remove(norm_order_betadata)
> gc()
          used (Mb) gc trigger   (Mb)   max used    (Mb)
Ncells  613789 32.8    5204299  278.0    8131716   434.3
Vcells 3295242 25.2 1121722249 8558.1 1752665689 13371.8
##20200423########################################
##test y1 -y8 province all selected top CPG　with Ｐ　＜　１ｅ－６
## try add genome PC as covar
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1'))  
library(data.table())
library(qqman)
workpath <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/"
setwd(workpath)
yname = "y1"
covnames = c("PC1", "PC2", "PC3", "PC4", "PC5", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age")
resultname <- "province.GABmale.conditionalPC.ewas1e-6.y1"
#ewas fun 
selectCPG <- matrix(NA, 20,4)
i = 1
	    phedata <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/usedCon2.GABmale985Methy.csv"
        savepath <- paste0(workpath, "y1PC5CellAgeresult_",i)
        ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
        #filter P < 1e-6 save
	    tmp <- paste0("awk -F, '$5<1e-6 {print $1,$7,$8,$5}' ",savepath, "/ewasresult_",yname,".csv"," > ",workpath,resultname,"_",i)
		system(tmp)
## 0CPG 
resultname <- "province.GABmale.conditionalPC.ewas5e-2.y1"
tmp <- paste0("awk -F, '$5<0.05 {print $1,$7,$8,$5}' ",savepath, "/ewasresult_",yname,".csv"," > ",workpath,resultname,"_",i)
		system(tmp)
#cg07315521,0.999597095762812
######################################
##select XinJiang top3 CPG
/LiuFanBin/R-3.6.1
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1'))  

workpath <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/"
setwd(workpath)
yname = "y1"
covnames = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age","cg15146515","cg03334052")
resultname <- "province.GABmale.conditional.ewas1e-6.y1"
i=3
phedata <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/province.GABmale.conditional.ewas1e-6.y1_2.phecondata_2"
savepath <- paste0(workpath, yname,"CellAgeresult_",i)
ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
#filter P < 1e-6 save
tmp <- paste0("awk -F, '$5<1e-6 {print $1,$7,$8,$5}' ",savepath, "/ewasresult_",yname,".csv"," > ",workpath,resultname,"_",i)
		system(tmp)
		# select top CPG!!!! NO.1 2  cg15146515 cg03334052 occur next by cg08597808
		# use code with i != 1
		covnames = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age","cg15146515","cg03334052")
		dataManhattan[which(dataManhattan$CPG %in% covnames),c(1:4)] <- NA
	    selectCPG[i,] <- as.matrix(dataManhattan[which(dataManhattan$P == min(na.omit(dataManhattan)$P)),c(1:4)])
######################################
##select XinJiang top4 CPG
workpath <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/"
setwd(workpath)
yname = "y1"
covnames = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age","cg15146515","cg03334052", "cg08597808")
resultname <- "province.GABmale.conditional.ewas1e-6.y1"
i=4
phedata <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/province.GABmale.conditional.ewas1e-6.y1_3.phecondata_3"
savepath <- paste0(workpath, yname,"CellAgeresult_",i)
ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
#Warning messages:
#1: glm.fit: fitted probabilities numerically 0 or 1 occurred
#filter P < 1e-6 save
###################################
###select XinJiang top5 CPG
covnames = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age","cg15146515","cg03334052", "cg08597808","cg00277334")
resultname <- "province.GABmale.conditional.ewas1e-6.y1"
i=5
phedata <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/province.GABmale.conditional.ewas1e-6.y1_4.phecondata_4"
savepath <- paste0(workpath, yname,"CellAgeresult_",i)
ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
#Warning messages:
# 1: glm.fit: fitted probabilities numerically 0 or 1 occurred

####################################
##select XinJiang top6 CPG
covnames = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "age","cg15146515","cg03334052", "cg08597808","cg00277334","cg24869232")
resultname <- "province.GABmale.conditional.ewas1e-6.y1"
i=6
phedata <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/province.GABmale.conditional.ewas1e-6.y1_5.phecondata_5"
savepath <- paste0(workpath, yname,"CellAgeresult_",i)
ewasFun(usenormbetadata = TRUE, phedata = phedata, selphe = yname, covariates = covnames, type = "logit", isbetadependent = FALSE, savepath = savepath, logitmaxit = 25, chunknum = 1, plotmq = TRUE)
# Warning messages:
#1: glm.fit: fitted probabilities numerically 0 or 1 occurred
# select top 6 "cg14514032"
# system('grep MemTotal /proc/meminfo')
# system('free -m')
############################################################################################
####20200426## try to predict province by chosen CPGS(y1-3 over)
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phedata1 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/province.GABmale.conditional.ewas1e-6.y1_7.phecondata_7", header = TRUE)
phedata2 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y2NeiMeng/province.GABmale.conditional.ewas1e-6.y2_4.phecondata_4", header = TRUE)
phedata3 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y3ShanXi/province.GABmale.conditional.ewas1e-6.y3_6.phecondata_6", header = TRUE)
phedata4 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y4HeNan/province.GABmale.conditional.ewas1e-6.y4_2.phecondata_2", header = TRUE)
phedata5 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y5ShanDong/province.GABmale.conditional.ewas1e-6.y5_1.phecondata_1", header = TRUE)
phedata6 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y6JiangXi/province.GABmale.conditional.ewas1e-6.y6_1.phecondata_1", header = TRUE)
phedata7 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y7SiChuan/province.GABmale.conditional.ewas1e-6.y7_1.phecondata_1", header = TRUE)
phedata8 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y8GuangXi/province.GABmale.conditional.ewas1e-6.y8_1.phecondata_1", header = TRUE)
data12 <- merge(phedata1, phedata2[,-c(2:25)], by = "Me.id")
data34 <- merge(phedata3[,-c(2:25)], phedata4[,-c(2:25)], by = "Me.id")
data56 <- merge(phedata5[,-c(2:25)], phedata6[,-c(2:25)], by = "Me.id")
data78 <- merge(phedata7[,-c(2:25)], phedata8[,-c(2:25)], by = "Me.id")
data1234 <- merge(data12, data34, by = "Me.id")
data5678 <- merge(data56, data78, by = "Me.id")
data18 <- merge(data1234, data5678, by = "Me.id")
write.csv(data18, file = "phecon23data1-8.csv", quote = FALSE, row.names = FALSE)
###20200427predict
phecondata <- read.csv("phecon23data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[c(3, 14:20)])
fo <- as.formula(paste0("province~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fit1 <- glm(fo, data = phecondata, family=binomial())
fitdata <- na.omit(phecondata)
library(nnet)
fit2 <- multinom(fo, data= na.oit(phecondata))
pre <- predict(fit2)
table(factor(na.omit(phecondata)$province), factor(pre))
library(caret)
per <- confusionMatrix(factor(pre),factor(fitdata$province))
library(pROC)
prob = predict(fit2, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
############################################
##20200507 predict 8 province 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phedata1 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y1XinJiang/province.GABmale.conditional.ewas1e-6.y1_7.phecondata_7", header = TRUE)
phedata2 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y2NeiMeng/province.GABmale.conditional.ewas1e-6.y2_4.phecondata_4", header = TRUE)
phedata3 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y3ShanXi/province.GABmale.conditional.ewas1e-6.y3_6.phecondata_6", header = TRUE)
phedata4 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y4HeNan/province.GABmale.conditional.ewas1e-6.y4_9.phecondata_9", header = TRUE)
phedata5 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y5ShanDong/province.GABmale.conditional.ewas1e-6.y5_6.phecondata_6", header = TRUE)
phedata6 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y6JiangXi/province.GABmale.conditional.ewas1e-6.y6_9.phecondata_9", header = TRUE)
phedata7 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y7SiChuan/province.GABmale.conditional.ewas1e-6.y7_9.phecondata_9", header = TRUE)
phedata8 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y8GuangXi/province.GABmale.conditional.ewas1e-6.y8_3.phecondata_3", header = TRUE)
data12 <- merge(phedata1, phedata2[,-c(2:25)], by = "Me.id")
data34 <- merge(phedata3[,-c(2:25)], phedata4[,-c(2:25)], by = "Me.id")
data56 <- merge(phedata5[,-c(2:25)], phedata6[,-c(2:25)], by = "Me.id")
data78 <- merge(phedata7[,-c(2:25)], phedata8[,-c(2:25)], by = "Me.id")
data1234 <- merge(data12, data34, by = "Me.id")
data5678 <- merge(data56, data78, by = "Me.id")
data18 <- merge(data1234, data5678, by = "Me.id")
write.csv(data18, file = "phecon53data1-8.csv", quote = FALSE, row.names = FALSE)
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[c(3, 14:20)])
fo <- as.formula(paste0("province~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
library(nnet)
fit2 <- multinom(fo, data= fitdata)
pre <- predict(fit2)
table(factor(na.omit(phecondata)$province), factor(pre))
library(caret)
per <- confusionMatrix(factor(pre),factor(fitdata$province))
library(pROC)
prob = predict(fit2, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc

#######################################################################################################
####20200510 select CpG annovar
###table_annovar.pl example/ex1.avinput humandb -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt
filepath=/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/test/y1CpGtoAnnova.txt
table_annovar=/thinker/storage/org/liufanGroup/liyi/software/annovar/table_annovar.pl
humandb=/thinker/storage/org/liufanGroup/liyi/software/annovar/humandb/
generef=/thinker/storage/org/liufanGroup/liyi/software/annovar/example/gene_fullxref.txt
cd /thinker/storage/org/liufanGroup/liyi/software/annovar
outname=/thinker/storage/org/liufanGroup/fanxiu/ancestry/GAB_1000G/Zang50HJD1000G/SNP.HJDZang.anno/ZangHJD2SNP_to_anno
perl $table_annovar  $filepath $humandb --outfile $outname --remove --buildver hg19 --protocol refGene,cytoBand,1000g2014oct_eur,1000g2014oct_afr,exac03,ljb26_all,clinvar_20140929,snp138 --operation gx,r,f,f,f,f,f,f --nastring . --csvout --polish  --xref $generef
## annovar fail because of no allele
##0511 used NCBI genome data viewer top cpg  gene
##try box/scatter plot 
## box-dotplot 散点和箱式图的合体
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
fitdata <- na.omit(phecondata)
y <- factor(ifelse(fitdata$y1==0,"ctrl","XinJiang"))
x <- fitdata$cg15146515
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg15146515")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
y <- factor(ifelse(fitdata$y2==0,"ctrl","NeiMeng"))
x <- fitdata$cg01083584
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg01083584")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
y <- factor(ifelse(fitdata$y3==0,"ctrl","ShanXi"))
x <- fitdata$cg20465143
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg20465143")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()

y <- factor(ifelse(fitdata$y4==0,"ctrl","HeNan"))
x <- fitdata$cg26517646
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg26517646")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
y <- factor(ifelse(fitdata$y5==0,"ctrl","ShanDong"))
x <- fitdata$cg00270311
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg00270311")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
y <- factor(ifelse(fitdata$y6==0,"ctrl","JiangXi"))
x <- fitdata$cg26104781
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg26104781")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
y <- factor(ifelse(fitdata$y7==0,"ctrl","SiChuan"))
x <- fitdata$cg25620050
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg25620050")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
y <- factor(ifelse(fitdata$y8==0,"ctrl","GuangXi"))
x <- fitdata$cg13580380
boxplot(x ~ y, border=4, col="light blue", boxwex=0.5, main = "cg13580380")
points(x ~ y,col=2,cex=0.5, pch=16)
dev.off()
##################################################################
awk -F, '{print $1,$7,$8,$5}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/y8GuangXi/y8CellAgeresult_4/ewasresult_y8.csv > /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/conManhattanPlot/usedManhattan.ewasresult_y8.csv

# /LiuFanBin/R-3.6.1 
##GABmale CpG cell7 age
## province 8 conditional GWAS
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1'))  
library(qqman)
library(data.table)
outpath <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/conManhattanPlot/"
setwd(outpath)
# read data
resultName <- "province.GABmale.conditional.ewas1e-6.y8"
yname <- "GuangXi"
dataManhattan <- fread("usedManhattan.ewasresult_y8.csv", header = TRUE)
names(dataManhattan) <- c("CPG", "CHR", "BP", "P")
conCpG <- read.table("y8.selectCpG.province.GABmale.conditional.ewas1e-6",header = TRUE)
# Manhattan plot
cpgnames <- conCpG$CPG
dataManhattan[which(dataManhattan$CPG %in% cpgnames),c(1:4)] <- NA
usedplot <- rbind(na.omit(dataManhattan), na.omit(conCpG))
SNPs <- as.character(na.omit(conCpG)$CPG)
names(usedplot) <- c("SNP", "CHR", "BP", "P")
tiff(paste0(yname,".ConditionalEWAS.anno.manhattan.tiff"), bg="white",height = 2400,width = 2600,res = 300)
manhattan(usedplot, col = c("red", "blue"),  highlight=SNPs, annotatePval=1e-6, annotateTop = FALSE, ylim=c(0, ceiling(-log10(min(usedplot$P)))), 
                    suggestiveline=FALSE, genomewideline = FALSE, main = paste0(yname,".Manhattan"))
abline(h=-log10(1e-6),col="red")
dev.off()
####################################################################################
##0512 leave1outAUC/auc
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[c(3, 14:20)])
fitdata <- na.omit(phecondata)
library(caret)

pre <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/test/province8.cpg53cell7age.pre", header = TRUE)
pre <- factor(pre$x)
table(factor(fitdata$province), pre)
per <- confusionMatrix(pre,factor(fitdata$province))
library(pROC)
prob = predict(fit2, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
###########################################################################
##20200524 leave1out province8~age +cpg53    (no cell 7)
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
# just cov age
covnames <- names(phecondata[3])
fitdata <- na.omit(phecondata)
library(caret)
fo <- as.formula(paste0("province~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
library(nnet)
fit2 <- multinom(fo, data= fitdata)
pre <- predict(fit2)
table(factor(na.omit(phecondata)$province), factor(pre))
library(caret)
per <- confusionMatrix(factor(pre),factor(fitdata$province))
library(pROC)
prob = predict(fit2, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
##waiting for pre and pro dataleave1out
pre <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/province8.cpg53age.pre", header = TRUE)
pre <- factor(pre$x)
table(factor(fitdata$province), pre)
per <- confusionMatrix(pre,factor(fitdata$province))
# result no cell 7  = +cell7
###########################################################################
##20200525 leave1out province8~age +cpg53 + population(9)    (no cell 7)
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
# test y1 ~ cpg7 +age +population
cpgnames <- names(phecondata[26:32])
# just cov age+cell7
covnames <- names(phecondata[c(3, 14:20)])
fitdata <- na.omit(phecondata)
fo <- as.formula(paste0("y1~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y2 ~ cpg4 +age +cell7
cpgnames <- names(phecondata[33:36])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y2~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y3 ~ cpg6 +age +cell7
cpgnames <- names(phecondata[37:42])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y3~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y4 ~ cpg9 +age +cell7
cpgnames <- names(phecondata[43:51])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y4~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y5 ~ cpg6 +age +cell7
cpgnames <- names(phecondata[52:57])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y5~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y6 ~ cpg9 +age +cell7
cpgnames <- names(phecondata[58:66])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y6~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y7 ~ cpg9 +age +cell7
cpgnames <- names(phecondata[67:75])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y7~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)
# test y8 ~ cpg3 +age +cell7
cpgnames <- names(phecondata[76:78])
# just cov age+cell7
covnames <- names(phecondata[c(3:4, 14:20)])
fo <- as.formula(paste0("y8~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
fit1 <- glm(fo, data = fitdata, family = binomial(link = "logit"))
summary(fit1)



###########################################################################
##20200525  plot pca plot by select 53cpg 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[c(3, 14:20)])
fitdata <- na.omit(phecondata[, c(5, 26:ncol(phecondata))]) #983
# set colours
PCphecon <- prcomp(as.matrix(fitdata[2:ncol(fitdata)]))
str(PCphecon)
# explained variance
PCphecon$sdev^2/sum(PCphecon$sdev^2)
# [1] 0.130201141 0.105172594
usedPC <- PCphecon$x
usedPC12 <- data.frame(usedPC[,c(1:2)])
usedPC12$pop <- factor(fitdata$province)
# 
require("RColorBrewer")
col <- colorRampPalette(c("red","orange", "yellow","forestgreen","grey","royalblue","purple","black"))(length(unique(usedPC12$pop)))[factor(usedPC12$pop)]
# generate PCA bi-plots
# project.pca <- eigenvec
# summary(project.pca)
# plot pc1-pc2
# par(mar=c(5,5,5,5), cex=2.0, cex.main=7, cex.axis=2.75, cex.lab=2.75, mfrow=c(1,2))
tiff("pc1.pc2.cpg53.province8.sample988.tiff", bg = "white", width = 3000, height = 2000, res = 240)
plot(usedPC12[,1], usedPC12[,2], type="n", main="pca.CpG53", adj=0.5, xlab="PC1(13.0%)", ylab="PC2(10.5%)", font=2, font.lab=2)
points(usedPC12[,1], usedPC12[,2], col=col, pch=20, cex=1.1)
legend("bottomright", bty="n", cex=1.1, title="", c("XinJiang","NeiMeng", "ShanXi", "HeNan", "ShanDong", "JiangXi", "SiChuan", "GuangXi"), fill=c("red","orange", "yellow","forestgreen","grey","royalblue","purple","black"))
dev.off()


## ## MeQTL(SNP ~select cpg53) predict province
cd /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result
head /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/cisSigRlt
head /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/111SNPtoProvince/select111SNP.list
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG
awk 'FNR==NR{a[$1]=$0;next}{ if(a[$1]) {print $0, a[$1]} }' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/select53CpG.list /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/cisSigRlt > meqtl.cismatch.cpg53.snp
awk 'FNR==NR{a[$1]=$0;next}{ if(a[$1]) {print $0, a[$1]} }' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/select53CpG.list /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/transSigRlt > meqtl.transmatch.cpg53.snp
## after select Cpg 360
awk '$4==cg15146515 {print $1, $2, $3, $4, $5}' /thinker/storage/org/liufanGroup/liuxx/MeQTL/whole_analyse_result/cisSigRlt > test1.cg15146515
#no meqtl CPG in LXX#
install.packages("MatrixEQTL")
##Warning: unable to access index for repository https://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib:
  cannot download all files
Warning messages:
1: In download.file(url, destfile = f, quiet = TRUE) :
  URL 'https://cran.r-project.org/CRAN_mirrors.csv': status was 'Couldn t resolve host name'
2: package ‘MatrixEQTL’ is not available (for R version 3.3.1)
### map cpg53 for MatrixEQTL but package fail !!!
plink --bfile /thinker/storage/org/liufanGroup/liufan/data/1000G/1000G --recode --out  
## try plink cpg53~SNP(GAB　inteersct3qc ) 
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53
#keep 988 sample genofile
plink --bfile /thinker/storage/org/liufanGroup/public_data/GABgenodata/Intersect_3qc_bim_SNP --keep /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/geo988list.usedmultilogis --make-bed --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988
# try cpg ~ genome cpg cg15146515  by plink 
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988  --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --pheno-name cg15146515 --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --covar-name age, B, NK, CD4T CD8T, Mono, Neutro, Eosino, PC1, PC2, PC3, PC4, PC5 --linear hide-covar --ci 0.95 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp
##but all th p (snp) = NA  select < 0.05 get zero
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53
awk '$12<0.05 {print $1, $2, $3, $4, $5, $6, $12}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.assoc.linear > cg15146515.meqtl.snp0.05
# try cov age cell 7 no pc1-5
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988  --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --pheno-name cg15146515 --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --covar-name age, B, NK, CD4T CD8T, Mono, Neutro, Eosino --linear hide-covar --ci 0.95 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covcell7age
awk '$12<0.05 {print $1, $2, $3, $4, $5, $6, $12}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covcell7age.assoc.linear > cg15146515.meqtl.snp.covcell7age.p0.05
# still no snp p< 0.05
# try cov just age 
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988  --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --pheno-name cg15146515 --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --covar-name age --linear hide-covar --ci 0.95 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covage
awk '$12<0.05 {print $1, $2, $3, $4, $5, $6, $12}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covage.assoc.linear > cg15146515.meqtl.snp.covage.p0.05
awk '$12 < 0.05/5983264 {print $1, $2, $3, $4, $5, $6, $12}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covage.assoc.linear > cg15146515.meqtl.snp.covage.pBonf0.05
# 47713 meQTL SNPs selected
# try cov age and pc1 +pc2
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988  --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --pheno-name cg15146515 --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --covar-name age, PC1, PC2 --linear hide-covar --ci 0.95 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covagepc12
awk '$12 < 0.05/5983264 {print $1, $2, $3, $4, $5, $6, $12}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg15146515.meqtl.snp.covagepc12.assoc.linear > cg15146515.meqtl.snp.covagepc12.pBonf0.05
# select 24 same chr snps  with cov ageand pc12
### so next use age +pc12 for other cpg52
######cg03334052
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988 --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --pheno-name cg03334052 --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --covar-name age, PC1, PC2 --linear hide-covar --ci 0.95 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg03334052.meqtl.snp.covagepc12
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988 --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --pheno-name cg15146515 --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis --covar-name age,B,NK,CD4T,CD8T,Mono,Neutro,Eosino --linear hide-covar --ci 0.95 --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53covagecell7/cg15146515_meqtl.snp.covagecell7
awk '$12 < 0.05/5983264 {print $1, $2, $3, $4, $5, $6, $12}' /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cg03334052.meqtl.snp.covagepc12.assoc.linear > cg03334052.meqtl.snp.covagepc12.pBonf0.05
##############################################################################
## try cpg53~age +cell7
genofile <- " --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988"
phenofile <- " --pheno /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis"
covarfile <- " --covar /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/plink.phecon53data1-8.usedmultilogis"
outpath <- "/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53covagecell7/"
outname <- "meqtl.covagecell7"
setwd(outpath)
phecondata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[c(3, 14:20)])
covname <- paste0(" --covar-name ",paste0(covnames,collapse = ","))
for (i in 1:length(cpgnames)){
     phename <- paste0(" --pheno-name ", cpgnames[i])
	 tmp <- paste0("plink",genofile,phenofile,phename,covarfile,covname," --linear hide-covar --ci 0.95"," --out ",outpath,cpgnames[i],"_",outname)
     system(tmp)
	 #filter P < 5e-8 save 
	 tmp <- paste0("awk '$12<0.05/5983264 {print $1, $2, $3, $4, $5, $6, $12}' ",outpath,cpgnames[i],"_",outname,".assoc.linear"," > ",outpath,cpgnames[i],"_",outname,"pBonf0.05")
	 system(tmp)	 
}


outname <- "meqtl.snp.covagecell7"
resultname <- "ancestry.intersect3qc.conditional.gwas.y2.result"

#############################################################
##20200601
##1.cpg_norm~age+3pc+3mePC+7cell ;2.cpg_norm_residuals~SNP
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull")
load("/thinker/storage/org/liufanGroup/liuxx/MeQTL/newQC/residuals_norm_order_betadata_ready.RData")
phecondata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
IID <- rownames(residuals_norm_order_betadata_ready)
selectCPG <- residuals_norm_order_betadata_ready[,which(colnames(residuals_norm_order_betadata_ready) %in% cpgnames)]
selectCPG <- data.frame(selectCPG)
selectCPG$IID <- IID
write.table(selectCPG, file="select53CPG.residuals_norm_order_betadata_ready", quote = FALSE, row.names = FALSE)
phedata <- phecondata[1:25]
residata <- merge(phedata, selectCPG, by.x ="Me.id", by.y="IID")
FID <- residata$Me.id
IID <- residata$Me.id
residata2 <- cbind(FID, residata)
write.table(residata2, file = "select53CPG.residata", quote = FALSE, row.names= FALSE)
#  nohup /LiuFanBin/R-3.6.1 < /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull/residcpg53.meqtl.covnull.R --no-save > /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull/residcpg53.meqtl.covnull.R.log &
#　[1] 59590
#################################################################
##20200602 meqtl CpG not in cpg53 %in% wrong
## sum(cpgnames1 %in% cpgnames)
#[1] 53
#> sum(cpgnames %in% cpgnames1)
#[1] 53
# select 20 snp meqtl ~ residata cpg53
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/genofile.GABgeo988 --extract /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull/snp20list.resid53cpg.meqtl --recodeA --out select20SNP.resid53cpg.meqtl

geno <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/plinkCpG53/cpg53residcovnull/select20SNP.resid53cpg.meqtl.raw", header = TRUE)
phecondata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/phecon53data1-8.csv", header  =TRUE)
geno <- geno[, -c(1,3,4,5,6)]
IID <- phecondata$Me.id
phecondata$IID <- IID
phecondata <- phecondata[, -1]
usedLogit <- merge(phecondata, geno, by = "IID")
str(usedLogit)

# rename colnames(problem from recodeA rs123 -> rs123_A)
str1 <- unlist(strsplit(names(usedLogit[79:ncol(usedLogit)]), "_"))
usedLogitSNP <- c()
i=1
while( i <= length(str1)){
usedLogitSNP <- c(usedLogitSNP,str1[i])
i= i+2
}
names(usedLogit)[79:ncol(usedLogit)] <- usedLogitSNP
## save usedLogit data
write.table(usedLogit, file = "province8.sample988.pheconmeQTLSNPdata.cpg53snp20", quote = FALSE, row.names = FALSE)

## predict 8 province by meqtl 20SNP
/LiuFanBin/R-3.6.1
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/residcpg53.meqtlsnp20")
data <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/meQTLs53CpG/residcpg53.meqtlsnp20/province8.sample988.pheconmeQTLSNPdata.cpg53snp20", header  =TRUE)
meqtlSNP <- names(data[79:ncol(data)])
library(nnet)

library(caret)
fo <- as.formula(paste0("province~", paste(meqtlSNP[1:length(meqtlSNP)], collapse = "+")))
fitdata <- na.omit(data)
fit <- multinom(fo, data= fitdata)
pre <- predict(fit)
per <- confusionMatrix(factor(pre),factor(fitdata$province))
library(pROC)
prob = predict(fit, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
###plot accumulate AUC plot for 8 province
## AUC plot::::::::::::::::::::::::::::::::::::::::::::::::: 
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
AUC <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/AUC.leave1out.province8.cpg53age", header = T)
numcpg <- nrow(AUC)
numy <- ncol(AUC)
number <- 1:numcpg
usedplot <- data.frame(number, AUC)
mydata<-usedplot
#color
color<-rainbow(numy)
rep <- rep(" ", times = numcpg)
#legend 
leg<-c("XinJiang", "NeiMeng",  "ShanXi", "HeNan", "ShanDong", "JiangXi", "SiChuan", "GuangXi")
# point shape
pch_vec <- rep(c(15:18), 2)
tiff("AUC increases as Number of CPGs increases.tiff",bg="white",height = 2400,width = 3000,res = 300)
plot(mydata[,1],mydata[,2],type = "o",pch=pch_vec[1],col=color[1],ylim = c(0.5,1),xaxt="n",xlab = "Number of CPGs",ylab = "AUC")
# axis side =1 = bottom
axis(side = 1,at = c(1:numcpg),labels = 1:numcpg, las = 0)
#text(x=c(1:56),y = 0.45, adj = 1, labels = 1:56,xpd = TRUE,xaxt="n")
for(i in 3:ncol(mydata)){
  lines(mydata[,1],mydata[,i],type = "o",pch=pch_vec[i-1],col=color[i-1],ylim = c(0,1),xaxt="n")
}
legend("bottomright",leg,pch = pch_vec,lty=rep(1,numy), border = TRUE,trace=TRUE,col = color,pt.cex = 1.7,lwd=1.5, bty = "o")
dev.off()
#################################################################
##20200619 try to anno 53 cpg
gm132-2
R
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/select53CpGtoAnno")
source("http://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
 .libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1'))

install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene")
##!!!fail host  name
#####################################################################
##20200620 meeeting Han(cpg vs snp)

setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition")
 phedata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/usedCon2.GABmale985Methy.csv",header = TRUE)
han <- phedata[which(phedata$population==1),]
write.csv(han, file = "usedCon2.GABhan495Methy.csv", quote = FALSE, row.names = FALSE)
#  3   4   5   6   7
# 150  97  50  96  98
## 20200621 select shanxi 6CpG
## next HeNan shandong jiangxi sichuan
####################################################################
##20200623 predict 8 province 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/HanMultilogis")
phedata3 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/y3ShanXi/province.GABmaleHan.conditional.ewas1e-6.y3_6.phecondata_6", header = TRUE)
phedata4 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/y4HeNan/province.GABmaleHan.conditional.ewas1e-6.y4_5.phecondata_5", header = TRUE)
phedata5 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/y5ShanDong/province.GABmaleHan.conditional.ewas1e-6.y5_2.phecondata_2", header = TRUE)
phedata6 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/y6JiangXi/province.GABmaleHan.conditional.ewas1e-6.y6_3.phecondata_3", header = TRUE)
phedata7 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/y7SiChuan/province.GABmaleHan.conditional.ewas1e-6.y7_4.phecondata_4", header = TRUE)
data34 <- merge(phedata3[,-c(2:25)], phedata4[,-c(2:25)], by = "Me.id")
data56 <- merge(phedata5[,-c(2:25)], phedata6[,-c(2:25)], by = "Me.id")
data3456 <- merge(data34, data56, by = "Me.id")
data34567 <- merge(phedata7, data3456, by = "Me.id")
write.csv(data34567, file = "Han.phecondata.csv", quote = FALSE, row.names = FALSE)
write.csv(na.omit(data34567), file = "Han483.phecondata.csv", quote = FALSE, row.names = FALSE)

setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/HanMultilogis")
phecondata <- read.csv("Han.phecondata.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[c(3, 14:20)])
fo <- as.formula(paste0("province~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
library(nnet)
fit2 <- multinom(fo, data= fitdata)
pre <- predict(fit2)
table(factor(na.omit(phecondata)$province), factor(pre))
library(caret)
per <- confusionMatrix(factor(pre),factor(fitdata$province))
library(pROC)
prob = predict(fit2, type='probs')
auc <- c(NA)
for (k in 1:5){
          #calculate 5 auc for each province 
		  roc_result <- roc(fitdata[, k+7], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
#####################################################################################
## 20200622 GABmale Han conditional GWAS
plink --bfile /thinker/storage/org/liufanGroup/public_data/GABgenodata/Intersect_3qc_bim_SNP --keep /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/list.QCusedCon2.GABmale483 --make-bed --out /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/genofile.GABmaleHan483
## phenofile y 0/1 -> 1/2/
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS")
data483 <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/province5.han.sample483.csv", header = TRUE)
write.table(data483, file = "province5.han.sample483", quote = F, row.names = F)
## after SNP select 
## construct multilogistic model
# recode 25 SNP to 0/1/2 -- unique
## R
cd /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis
plink --bfile /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/genofile.GABmaleHan483 --extract /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis/select25SNP.list --recodeA --out select25SNP.GABmaleHan483

setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis")
geno <- read.table("select25SNP.GABmaleHan483.raw", header = TRUE)
province <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/province5.han.sample483", header = TRUE)
geno <- geno[, -c(1,3,4,5,6)]
province <- province[, -1]
usedLogit <- merge(province, geno, by = "IID")
str(usedLogit)
# replace X to " ", . to ":"(problem from read.table .raw [8: -> X8.])
names(usedLogit) <- gsub("X", "", names(usedLogit))
names(usedLogit) <- gsub("\\.", ":", names(usedLogit))
# set variable class
usedLogit$y1 <- factor(usedLogit$y1)
usedLogit$y2 <- factor(usedLogit$y2)
usedLogit$y3 <- factor(usedLogit$y3)
usedLogit$y4 <- factor(usedLogit$y4)
usedLogit$y5 <- factor(usedLogit$y5)
usedLogit$y6 <- factor(usedLogit$y6)
usedLogit$y7 <- factor(usedLogit$y7)
usedLogit$y8 <- factor(usedLogit$y8)
usedLogit$province <- factor(usedLogit$province)
# rename colnames(problem from recodeA rs123 -> rs123_A)
str1 <- unlist(strsplit(names(usedLogit[14:ncol(usedLogit)]), "_"))
usedLogitSNP <- c()
i=1
while( i <= length(str1)){
usedLogitSNP <- c(usedLogitSNP,str1[i])
i= i+2
}
names(usedLogit)[14:ncol(usedLogit)] <- usedLogitSNP
## save usedLogit data
save(usedLogit, file = "provinceHan.sample483.usedLogit.RData")
write.table(usedLogit, file = "provinceHan.sample483.usedLogit", quote = FALSE, row.names = FALSE)
## run multilogistic predict
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis")
library(dplyr)
library(caret)
library("nnet")
load("provinceHan.sample483.usedLogit.RData")
fitdata <- na.omit(usedLogit[,-c(1:4)]) 
fitdata <- na.omit(fitdata[,-c(2:9)]) 
logis1 <- multinom(province~., fitdata)
per_logis1 <-  confusionMatrix(factor(predict(logis1)),fitdata$province)
per_logis1
## AUC
prob <-  predict(logis1, newdata = fitdata[,-1], type='probs')
auc <- NA
usedLogit <- na.omit(usedLogit[, -c(1:4)]) 
for (i in 1:5){
    roc_result <- roc(usedLogit[, i+3], prob[,i], auc=T, quiet = TRUE)
    auc[i] <- roc_result$auc
}
auc
######################################
##20200626 leave1outAUC conditionGWAS(SNP25) VS conditionalEWAS(CpG20)
#!! SNP name 1:1213 -> X1.1213 invlid name for model
##cpg 
#################################################################
#20200701:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#add Han CpG (vs SNP leave1out result) 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/HanMultilogis")
library(caret)
phecondata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/HanMultilogis/Han483.phecondata.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
# cov age
covnames <- names(phecondata[3]) 
fitdata <- na.omit(phecondata)
yname <- as.character(paste("y", c(3:7), sep = ""))
pre <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConEWAS/HanMultilogis/accumulateCpGtoProvince5Han.pre", header = TRUE)
# !!prevince3-7 predict to 1-5  so pre$x+2
per_logis1 <-  confusionMatrix(factor(pre$x+2),factor(fitdata$province))
per_logis1

#add Han CpG (vs SNP leave1out result) 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis")
library(caret)
usedLogit <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis/provinceHan.sample483.usedLogit.csv", header = TRUE)
fitdata <- na.omit(usedLogit[,-c(1:4)]) #475*34
snpname <- colnames(fitdata)[10:ncol(fitdata)]
yname <- as.character(paste("y", c(3:7), sep = ""))
pre <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/HanCondition/hanConGWAS/HanMultilogis/accumulateSNPtoProvince5Han.pre", header = TRUE)
# !!prevince3-7 predict to 1-5  so pre$x+2
per_logis1 <-  confusionMatrix(factor(pre$x+2),factor(fitdata$province))
per_logis1
######################################################
##test CpG 53 pre geo AREA 7
phecondata <- read.csv("phecon53geo7.csv", header  =TRUE)
# add col geo7 so cpgnames from col 26+1(geo7)-1(y8)=26
cpgnames <- names(phecondata[26:ncol(phecondata)])
# coa age 965*79
covnames <- names(phecondata[3]) 
fo <- as.formula(paste0("geo7~", paste(cpgnames[1:length(cpgnames)], collapse = "+"),"+", paste(covnames[1:length(covnames)], collapse = "+")))
fitdata <- na.omit(phecondata)
library(nnet)
fit2 <- multinom(fo, data= fitdata)
pre <- predict(fit2)
library(caret)
per <- confusionMatrix(factor(pre),factor(fitdata$geo7))
library(pROC)
prob = predict(fit2, type='probs')
auc <- c(NA)
for (k in 1:7){
          #calculate 7 auc for each province 
		  roc_result <- roc(fitdata[, k+6], prob[,k], auc=T, quiet = TRUE) ##y1=factor(y1)=y1-1
          auc[k] <- roc_result$auc
          } 
auc
#20200702:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(caret)
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/geo7pre")
pre <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/geo7pre/geoarea7.cpg53age.pre", header = TRUE)
fitdata <- na.omit(phecondata)
per_logis1 <-  confusionMatrix(factor(pre$x),factor(fitdata$geo7))
per_logis1
######################################################################
##20200702 infer_G software::::::::::::::::::::::::::::::::::::::::::
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/inferG_software")
phecondata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
covnames <- names(phecondata[3]) 
fitdata <- na.omit(phecondata)#965*78
age <- fitdata$age
y <- fitdata$province
inferG <- data.frame(y, fitdata[,c(26:ncol(fitdata))])
inferG$y <- factor(inferG$y)
save(inferG, file = "inferG.RData")
write.csv(inferG, file = "inferG.csv", quote = FALSE, row.names= FALSE)
##deleat cov age 
y ~cpg53 [no cov]

## ########################################
##20200703 share GAB Methylation data 998 after QC to JiangLi and XuJiChen　
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/toGABMethyData")
load("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/betadata.RData")
ls()
#[1] "betadata"
dim(betadata) 
# [1]    988 723730
betadata1 <- betadata[c(1:500),]
save(betadata1, file = "betadata1.RData")
rownames(betadata1)
betadata2 <- betadata[c(501:nrow(betadata)),]
rownames(betadata2)
save(betadata2, file = "betadata2.RData")
####################################################################
##20200714 add AUC for cov Null (-age -cell)
library(caret)
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
phecondata <- read.csv("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/phecon53data1-8.csv", header  =TRUE)
pre <- read.table("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/province8.cpg53covNull.pre", header = TRUE)
fitdata <- na.omit(phecondata)
per_logis1 <-  confusionMatrix(factor(pre$x),factor(fitdata$province))
per_logis1
######################################################################
##20200726 ANOVA of cpg53 and mean 
cpgdata <- fitdata[,c(5,26:ncol(fitdata))]
cpgnames <- names(fitdata)[26:ncol(fitdata)]
# cpgdata <- as.matrix(cpgdata) group _by must be data.frame
# cpgdata2 <- group_by(cpgdata,province)
# mean1 <- summarise(cpgdata2, avg_cg15146515 = mean(cg15146515)) 
# cpg <- cpgnames[i]
# avg <- summarise(cpgdata2, avg = mean()) 
##!! group_by + mean -> for loop fail turn to fun(aggregate)
avgCpG <- matrix(NA, length(cpgnames), 8)
for (i in 1:53){
  avgCpG[i,] <- aggregate(cpgdata[,i+1], by = list(cpgdata$province), FUN = mean)[,2]
}
rownames(avgCpG) <- cpgnames
write.table(avgCpG, file = "avg.cpg53.province8", quote = FALSE)
########################################################################
##20200801select 53CpG CHR BP 200kb
ordata <- data[order(data$CHR,data$BP),]
C500kb <- NA
# chr ==15 just 1 SNP
for (i in 1:22){
      numCSNP= sum(ordata$CHR==i)
	  tmp <- 0
	  if(numCSNP>2){
	        for (j in 2:numCSNP){
	               if((ordata[j,5] - ordata[j-1,5]) < 200000) #500kb
         	             {tmp <- tmp+1}
	  }
	  }
	  C500kb[i] <- tmp
	  ordata <- ordata[-c(1:numCSNP),] 
}