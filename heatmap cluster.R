## Heatmap clust
library(gplots)
data(mtcars)
x <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start=0, end=.3)
cc <- rainbow(ncol(x), start=0, end=.3)
    Cluster_Method<-c( "ward.D", "single", "complete", "average", "mcquitty", "median", "centroid")
    #R语言里面自带的hclust函数共有7种聚类方法
    for (i in 1:length(Cluster_Method)){
     #make a function to extract the cluster method
     myclust<-function(x){
     hclust(x,method=Cluster_Method[i])
     }
     #make heatmap by jpeg
     jpeg(filename=paste(Cluster_Method[i],'.jpg'),width=1024,height=728)
     heatmap.2(x,
     trace='none',
     hclustfun=myclust,
     labRow=NA,
     xlab='CellLines',
     ylab='Probes',
     main=Cluster_Method[i],
     col=greenred(64))
     dev.off()
    }
heatmap.2(x)
heatmap.2(x, dendrogram="none")
    heatmap.2(x, dendrogram="row") # 只显示行向量的聚类情况
    heatmap.2(x, dendrogram="col") #只显示列向量的聚类情况
    heatmap.2(x, srtCol=NULL)
    heatmap.2(x, srtCol=0, adjCol = c(0.5,1) )
    heatmap.2(x, srtCol=45, adjCol = c(1,1) )
    heatmap.2(x, srtCol=135, adjCol = c(1,0) )
    heatmap.2(x, srtCol=180, adjCol = c(0.5,0) )
    heatmap.2(x, srtCol=225, adjCol = c(0,0) ) ## not very useful
    heatmap.2(x, srtCol=270, adjCol = c(0,0.5) )
    heatmap.2(x, srtCol=315, adjCol = c(0,1) )
    heatmap.2(x, srtCol=360, adjCol = c(0.5,1) )
    heatmap.2(x, srtRow=45, adjRow=c(0, 1) )
    heatmap.2(x, srtRow=45, adjRow=c(0, 1), srtCol=45, adjCol=c(1,1) )
    heatmap.2(x, srtRow=45, adjRow=c(0, 1), srtCol=270, adjCol=c(0,0.5) )
    ## Show effect of offsetRow/offsetCol (only works when srtRow/srtCol is 
    ## not also present) heatmap.2(x, offsetRow=0, offsetCol=0)
    heatmap.2(x, offsetRow=1, offsetCol=1) 
    heatmap.2(x, offsetRow=2, offsetCol=2) 
    heatmap.2(x, offsetRow=-1, offsetCol=-1) 
    heatmap.2(x, srtRow=0, srtCol=90, offsetRow=0, offsetCol=0)
    heatmap.2(x, srtRow=0, srtCol=90, offsetRow=1, offsetCol=1)
    heatmap.2(x, srtRow=0, srtCol=90, offsetRow=2, offsetCol=2)
    heatmap.2(x, srtRow=0, srtCol=90, offsetRow=-1, offsetCol=-1)
    ## Show effect of z-score scaling within columns, blue-red color scale 
    hv <- heatmap.2(x, col=bluered, scale="column", tracecol="#303030")
    > names(hv) # 可以看到hv对象里面有很多子对象
    > "rowInd" "colInd" "call" "colMeans" "colSDs" "carpet" "rowDendrogram" "colDendrogram" "breaks" "col" "vline" "colorTable" ## Show the mapping of z-score values to color bins hv$colorTable
    ## Extract the range associated with white 我们得到了热图的颜色的数值映射矩阵，接下来就可以进行一系列的操作~！！！
    hv$colorTable[hv$colorTable[,"color"]=="#FFFFFF",]
    whiteBin <- unlist(hv$colorTable[hv$colorTable[,"color"]=="#FFFFFF",1:2]) 
    rbind(whiteBin[1] * hv$colSDs + hv$colMeans, whiteBin[2] * hv$colSDs + hv$colMeans )



