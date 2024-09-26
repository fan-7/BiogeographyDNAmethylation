########################################################
#'@Date 2019-09-23
#'@Author Jingxiaoxi and Liuxinxuan
#'@Description This is a function to do epigenome-wide association study (EWAS) on the 850K methylation data(all males) with 3 default files and 1 user supplied file, generates a .csv file contains EWAS result, and four plots (manhattan and qq). It has 6 parameters.
#'@Environment R; R packages ("plyr","qqman"); 3 default files in "/thinker/storage/org/liufanGroup/liuxx/DNAm_data/": betadata.RData (contains the 850K methylation beta values data, all males), genoinfo.csv (contains genomic information of the 850K probes) and defcov.csv (contains default covariates); this script applies to the sever:192.168.118.92.
#'@Param phedata: A csv file, giving phenotype data, contains a column named "Me.id" (individual id, factor or character, consistent with "Me.id" in defcov.csv) and interested phenotype datas (it can contain covariates, and any individual containing missing value in selected phenotype or covariates will be excluded in following analysis. If any phenotype in your phedata have the same name with any one in defcov.csv, the former will cover the latter's data in following analysis).
#'@Param selphe: String, giving phenotype's name selected to analyse. Must contained by phedata.
#'@Param covariates: Character vector or NULL, indicating which phenotypes as covariates. Default values are all phenotypes in defcov.csv, ie predicted blood cell fraction ("B", "NK", "CD4T", "CD8T", "Mono", "Neutro"), and the top five genomic principle components ("PC1", "PC2", "PC3", "PC4", "PC5"). You can also add other phenotypes as covariates by add them to your phedata and change the value of the covariates parameter.
#'@Param type: analytical method, "linear" (linear regression) or "logit" (logit regression). Default is "linear".
#'@Param savepath: String, giving path to save output files. Default is ".".
#'@Param logitmaxit: Integer, giving the maximal number of IWLS iterations when using logit regression. Default is 25.
#'@Return A piece of message like "page.mean ewas analytic (linear) results are saved in ./ewasFun_output".
#'@Export ewasresult_"selphe".csv, columns: probid (methylation probe id), Beta (beta value in regression), SE (standard error), STAT (statistic), P (p value), R.squared, FDR.P (p value adjusted by Benjamini-Hochberg procedure to control FDR), SNP, CHR, BP (base position), N (number of individuals analyzed) and gc.P (P value after genomic control).
#'@Export 4 plots, manhattan and qq plots before and after genomic control.
#'@FunctionsCalled R functions: merge, plyr::join, plyr::arrange, sapply, lm, glm, subset, tiff, suppressPackageStartupMessages, library, qqman::manhattan, qqman::qq, assign, detach.
#'@Example takes about 78 minutes
#'ewasFun(phedata = "./test.csv", selphe = "page.mean", covariates = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "PC1", "PC2", "PC3", "PC4", "PC5", "age"), type = "linear", savepath = "./ewasFun_output")
#'@Notice the primary loading file, "betadata.RData", is about 5GB size, so prepare enough memory before using the function.
#'@Notice if the phenotypes in your phedata have the same names with the ones in defcov.csv, the data in your phedata will be reserved for subsequent analysis, and the data in defcov.csv will be excluded.
#'@Notice by default, the selected phenotype is specified as the dependent variable and DNA methylation ¦Â value at each CpG site as the predictor of interest. You can modify the script in line 46 to change that.
########################################################
ewasFun <- function(phedata, selphe, covariates = c("B","NK","CD4T","CD8T","Mono","Neutro","Eosino","PC1","PC2","PC3","PC4","PC5","age","people","BMI"), type = "linear", savepath = ".", logitmaxit = 25)
{
    if(length(selphe) != 1)
    {
        stop("selphe must be a single string")
    }
    if(!dir.exists(savepath))
    {
        dir.create(file.path(savepath), recursive = TRUE)
    }
    suppressPackageStartupMessages(library(qqman))
    suppressPackageStartupMessages(library(plyr))
    phedata <- read.csv(file = phedata, as.is = TRUE)
    defcov <- read.csv("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/defcov.csv", as.is = TRUE) #default covariates
    samephe <- names(defcov) %in% names(phedata)
    samephe[1] <- FALSE
    defcov <- subset(defcov, select = !samephe)
    phedata <- merge(defcov, phedata, by = "Me.id") #merge user supplied phedata
    if(!all(c(selphe, covariates) %in% names(phedata)))
    {
        stop("selphe and covariates must be contained in the defcov.csv or phedata")
    }
    genoinfo <- read.csv("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/genoinfo.csv", as.is = TRUE) #genomic information
    load("/thinker/storage/org/liufanGroup/liuxx/DNAm_data/betadata.RData") #850K beta values data
    phedata <- na.omit(subset(phedata, select = c("Me.id", selphe, covariates)))
    betadata <- betadata[match(phedata$Me.id, rownames(betadata), nomatch = 0),] #select and arrange betadata
    phedata <- join(data.frame(Me.id = rownames(betadata)), phedata, by = "Me.id") #select and arrange phedata in the order of id in betadata
    phedata <- subset(phedata, select = c(selphe, covariates))
    fo <- formula(paste(selphe,"~ .")) #you can use fo <- formula("x ~ .") to specify each beta value as the dependent variable and the phenotype as the predictor of interest.
    
    #linear regression
    if(type == "linear")
    {
        linearfun <- function(x)
        {
            tmp <- cbind(x, phedata)
            smodel <- summary(lm(fo, data = tmp))
            c(smodel$coef[2,], R.squared = smodel$r.squared)
        }
        ewasresult <- data.frame(t(sapply(betadata, linearfun)))
    }
    #logit regression
    if(type == "logit")
    {
        logitfun <- function(x)
        {
            tmp <- cbind(x, phedata)
            summary(glm(fo, data = tmp, family = binomial(link = "logit"), control = list(maxit = logitmaxit)))$coef[2,]
        }
        ewasresult <- data.frame(t(sapply(betadata, logitfun)))
    }

    #rename and arrange the ewasresult
    names(ewasresult)[1:4] <- c("Beta", "SE", "STAT", "P")
    ewasresult$probid <- rownames(ewasresult)
    rownames(ewasresult) <- NULL
    ewasresult$FDR.P <- p.adjust(ewasresult$P, method = "fdr")
    ewasresult <- merge(ewasresult, genoinfo, by = "probid", all.x = TRUE)
    ewasresult$N <- nrow(phedata)
    ewasresult <- arrange(ewasresult, P)

    #gc genomic_control--Correction for Population Stratification
    ewasresult$gc.P <- NA
    n = nrow(betadata)
    tSTAT <- ewasresult$STAT[!is.na(ewasresult$P)]
    tN <- ewasresult$N[!is.na(ewasresult$P)]
    expected_stat <- rt(length(tSTAT), n)
    genomic_control <- median(abs(tSTAT)) / median(abs(expected_stat))
    tSTAT <- tSTAT / genomic_control
    ewasresult$gc.P[seq_along(tSTAT)][tSTAT < 0] <- pt(tSTAT[tSTAT < 0], tN[tSTAT < 0] - 1) * 2
    ewasresult$gc.P[seq_along(tSTAT)][tSTAT > 0] <-(1 - pt(tSTAT[tSTAT > 0], tN[tSTAT > 0] - 1)) * 2

    #save results
    #EWAS result
    write.csv(ewasresult, file = paste0(savepath, "/ewasresult_", selphe, ".csv"), quote = FALSE, row.names = FALSE)
    #manhattan and qq plots before gc
    tiff(filename = paste0(savepath, "/manhattan_plot_", selphe, ".tif"))
    manhattan(x = ewasresult, col = c("aquamarine4", "tomato"), suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08), main = paste0("manhattan_plot_", selphe))
    dev.off()
    tiff(filename = paste0(savepath, "/qq_plot_", selphe, ".tif"))
    qq(ewasresult$P, main = paste0("qq_plot_", selphe))
    dev.off()
    #manhattan and qq plots after gc
    tiff(filename = paste0(savepath, "/manhattan_gc_plot_", selphe, ".tif"))
    manhattan(x = cbind(subset(ewasresult, select = c("SNP", "BP", "CHR")), P = ewasresult$gc.P), col = c("aquamarine4", "tomato"), suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08), main = paste0("manhattan_gc_plot_", selphe))
    dev.off()
    tiff(filename = paste0(savepath, "/qq_gc_plot_", selphe, ".tif"))
    qq(ewasresult$gc.P, main = paste0("qq_gc_plot_", selphe))
    dev.off()
    detach("package:qqman")
    message(selphe, " ewas analytic (", type, ") results are saved in ", savepath)
}
#Example
ewasFun(phedata = "/thinker/storage/org/liufanGroup/qianyu/GAB2018/GH/sigSNP_cpg/sigSNP_cpg_GWAS.csv", selphe = "rs1020707975", covariates = c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino", "PC1", "PC2", "PC3", "PC4", "PC5", "age","people","BMI"), type = "linear", savepath = "/thinker/storage/org/liufanGroup/qianyu/GAB2018/GH/sigSNP_cpg")