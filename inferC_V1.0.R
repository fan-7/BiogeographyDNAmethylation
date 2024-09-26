# Source code for Chinese language family prediction tool 
# version: 1.0

geto <- function(df,newdata,tr=tr){
	# df: internal data consisting of all samples
	# newdata: if missing, infile is read for prediction
	# tr: threshold for AUC loss, predicted probabilities will not be reported if AUC loss > tr
	require(nnet)
	require(pROC)
	require(caret)
	snpinfo <- names(df)[-1]
	dftrain <- df[!is.na(df$y),names(df) %in% snpinfo | names(df) %in% "y"]
	fo <- as.formula(paste("y~",paste(snpinfo,collapse="+")))

	getauc <- function(fo,data,ycol = 1){
		fit <- multinom(fo=fo,data=data)
		yp <- predict(fit,type="probs")
		if(is.null(dim(yp))){yp<-as.matrix(yp)}
		ny <- ncol(yp)
		auc <- rep(NA,ny)
		for (j in 1:ny){
		    yt <- ifelse(data[,ycol]==j,1,0) 
			roc_result <- roc(yt, yp[,j],auc=T, quiet = TRUE)
            auc[j] <- roc_result$auc  # AUC
		}
		return (auc)
	}
	auc <- getauc(fo,dftrain) 
	
	getsenspe <- function(fo, data, ycol = 1){
	    fit <- multinom(fo=fo,data=data)
		yp <- predict(fit)
		if(is.null(length(yp))){yp<-as.matrix(yp)}
		per_fit <-  confusionMatrix(yp,dftrain$y)
		senspe <- matrix(NA, 2, nrow(per_fit$byClass))
		#Sensitivity
		senspe[1,] <- per_fit$byClass[,1] 
		#Specificity
        senspe[2,] <- per_fit$byClass[,2] 
		return (senspe)
	}
    senspe <- getsenspe(fo,dftrain)	 		

	# prediction result
	o <- list() 
	o[[1]] <- matrix(NA,dim(newdata)[1],26)  
	colnames(o[[1]])<-c(
		"CNH",
		"CSH",
		"TUR",
		"VIE",
		"MON",
		"HMO",
		"TIB",
		"KAM",
		"Full_AUC_CNH",
		"Full_AUC_CSH",
		"Full_AUC_TUR",
		"Full_AUC_VIE",
		"Full_AUC_MON",
		"Full_AUC_HMO",
		"Full_AUC_TIB",
		"Full_AUC_KAM",
		"Numb_missingSNPs",
		"Name_missingSNPs",
		"AUC_Loss_CNH",
		"AUC_Loss_CSH",
		"AUC_Loss_TUR",
		"AUC_Loss_VIE",
		"AUC_Loss_MON",
		"AUC_Loss_HMO",
		"AUC_Loss_TIB",
		"AUC_Loss_KAM"
	)
	o[[1]] <- as.data.frame(o[[1]])

	for (i in 1:nrow(newdata)){
		ix <- snpinfo %in% names(newdata)[!is.na(newdata[i,])]
		if(sum(ix)>0){
			fo <- as.formula(paste("y~",paste(snpinfo[ix],collapse="+")))
			fit <- multinom(fo,data=df) 
			o[[1]][i,1:8]<-predict(fit,newdata=newdata[i,],type="probs")
			o[[1]][i,9:16]<-auc
			o[[1]][i,17]<-sum(!ix)
			o[[1]][i,18]<-paste(snpinfo[!ix],collapse="/")
			aucnow <- getauc(fo,df)
			o[[1]][i,19:26]<-auc-aucnow
		}	
	}
	if(i==1){
		    if(any(o[[1]][1,19:26] > tr)){o[[1]][1,1:8]<-NA} 
	}else{
		ix <- apply(o[[1]][,19:26]>tr,1,function(x)any(x,na.rm=T))
		o[[1]][ix,1:8] <- NA
	}

	# output
	out <- newdata[,!names(newdata) %in% snpinfo,drop=F]
	out <- cbind(out,o[[1]])
	return(out)
}
# load internal data
load("inferC.RData")
# input newdata
input <- read.csv("input.csv", header = TRUE)
# run 
result <- geto(df=inferC,newdata = input, tr=0.5)
# output predict result
write.csv(result,"Result.csv",row.names=F,quote=F)

