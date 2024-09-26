## get accumulate AUC with cpg increase 
##2020713 multinom logistic regression 
## GAB972 Nna prov8  cov Null
## test leave one out 
#nohup /LiuFanBin/R-3.6.1 < /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/province8.CpG53.covNULL.leave1outAUC.R --no-save > /thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro/province8.CpG53.covNULL.leave1outAUC.R.log &
.libPaths(c(.libPaths(), '/thinker/storage/org/liufanGroup/public_software/R/3.6.1')) 
setwd("/thinker/storage/org/liufanGroup/fanxiu/geography/GAB1000Methy/condition2/multilogis8pro")
library(dplyr)
library("nnet")
library(pROC)
phecondata <- read.csv("phecon53data1-8.csv", header  =TRUE)
cpgnames <- names(phecondata[26:ncol(phecondata)])
# coa null
fitdata <- na.omit(phecondata)
yname <- as.character(paste("y", c(1:8), sep = ""))
auc <- matrix(c(NA), length(cpgnames), length(yname))
pre <- c(NA)
for (i in 1:length(cpgnames)){
     fo <- as.formula(paste0("province~", paste(cpgnames[1:i], collapse = "+")))
     prob <- matrix(NA, nrow(fitdata), length(yname))
     for(j in 1:nrow(fitdata)){
          #leave one for validation others for training
          validation<-fitdata[j,]
          training<- fitdata[-j,]
		  # fit model 
          logis1 <- multinom(fo, data = training, trace = FALSE)
		  # predict 8 probility 
		  prob[j,] = predict(logis1, newdata = validation, type='probs')
		  pre[j] <- predict(logis1, newdata = validation)
		  write.table(prob, file= "province8.cpg53covNull.pro", row.names = FALSE, quote = FALSE)
		  write.table(pre, file= "province8.cpg53covNull.pre", row.names = FALSE, quote = FALSE)
		  }
     for (k in 1:length(yname)){
          #calculate 8 auc for each language family class 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE)
          auc[i,k] <- roc_result$auc
          } 
write.table(auc, file = "AUC.leave1out.province8.cpg53.covNull", row.name = F, quote = F)
    }
write.table(auc, file = "AUC.leave1out.province8.cpg53.covNull", row.name = F, quote = F)