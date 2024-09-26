#################################################################
##China province inference by CPG
##accumulate AUC plot 
library(pROC)
prob = predict(fit, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each province 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) 
          auc[k] <- roc_result$auc
          } 
write.table(auc, file = "AUC.province8.cpg53age", row.name = F, quote = F)
# plot accumulate AUC plot for 8 province
AUC <- read.table("AUC.province8.cpg53age", header = T)
numcpg <- nrow(AUC)
numy <- ncol(AUC)
number <- 1:numcpg
usedplot <- data.frame(number, AUC)
mydata<-usedplot
# color
color<-rainbow(numy)
rep <- rep(" ", times = numcpg)
# legend 
leg<-c("XinJiang", "NeiMeng",  "ShanXi", "HeNan", "ShanDong", "JiangXi", "SiChuan", "GuangXi")
# point shape
pch_vec <- rep(c(15:18), 2)
tiff("AUC increases as Number of CPGs increases.tiff",bg="white",height = 2400,width = 3000,res = 300)
plot(mydata[,1],mydata[,2],type = "o",pch=pch_vec[1],col=color[1],ylim = c(0.5,1),xaxt="n",xlab = "Number of CPGs",ylab = "AUC")
# axis bottom
axis(side = 1,at = c(1:numcpg),labels = 1:numcpg, las = 0)
for(i in 3:ncol(mydata)){
  lines(mydata[,1],mydata[,i],type = "o",pch=pch_vec[i-1],col=color[i-1],ylim = c(0,1),xaxt="n")
}
legend("bottomright",leg,pch = pch_vec,lty=rep(1,numy), border = TRUE,trace=TRUE,col = color,pt.cex = 1.7,lwd=1.5, bty = "o")
dev.off()
