#################################################################
##China language family inference by SNP
##accumulate AUC plot 
library(pROC)
prob = predict(fit, type='probs')
auc <- c(NA)
for (k in 1:8){
          #calculate 8 auc for each language family 
		  roc_result <- roc(fitdata[, k+5], prob[,k], auc=T, quiet = TRUE) 
          auc[k] <- roc_result$auc
          } 
write.table(auc, file = "AUC.language8.111SNP", row.name = F, quote = F)
# plot accumulate AUC plot for 8 language family
AUC <- read.table("AUC.language8.111SNP", header = T)
numsnp <- nrow(AUC)
numy <- ncol(AUC)
number <- 1:numsnp
usedplot <- data.frame(number, AUC)
mydata<-usedplot
#color
color<-rainbow(numy)
rep <- rep(" ", times = numsnp)
#legend 
leg<-c("NorthHan", "SouthHan",  "Turkic", "Vietic", "Mongolia", "Hmong-Mien", "Tibeto-Burman", "Kam-Tai")
# point shape
pch_vec <- rep(c(15:18), 2)
tiff("AUC increases as Number of SNPs increases.tiff",bg="white",height = 2400,width = 3000,res = 300)
plot(mydata[,1],mydata[,2],type = "o",pch=pch_vec[1],col=color[1],ylim = c(0.5,1),xaxt="n",xlab = "Number of SNPs",ylab = "AUC")
# axis side =1 = bottom
axis(side = 1,at = c(1:numsnp),labels = 1:numsnp, las = 0)
#text(x=c(1:56),y = 0.45, adj = 1, labels = 1:56,xpd = TRUE,xaxt="n")
for(i in 3:ncol(mydata)){
  lines(mydata[,1],mydata[,i],type = "o",pch=pch_vec[i-1],col=color[i-1],ylim = c(0,1),xaxt="n")
}
legend("bottomright",leg,pch = pch_vec,lty=rep(1,numy), border = TRUE,trace=TRUE,col = color,pt.cex = 1.7,lwd=1.5, bty = "o")
dev.off()