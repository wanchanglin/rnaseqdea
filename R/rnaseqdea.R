#' Functions for RNA-Seq differential expression analysis
#' wl-23-06-2014, Mon: Commence
#' wl-29-06-2019, Sat: Tidy up

library(plyr)
library(reshape)
library(latticeExtra)
library(WriteXLS)
library(DESeq2)
library(DESeq)
library(edgeR)
library(NBPSeq)
library(EBSeq)
library(NOISeq)
library(samr)

#' ========================================================================
#' wll-14-08-2014: NGS count data analysis.
#' wll-15-08-2014: Return all plots individually.
#' wll-17-08-2014: Re-write boxplot
#' ========================================================================
ngs <- function(data, cls, com, method="stats.TSPM", norm.method="TMM"){
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  method <- 
    match.arg(method, 
              c("stats.DESeq","stats.DESeq2","stats.edgeR","stats.TSPM",
                "stats.NBPSeq","stats.EBSeq","stats.Wilcox","stats.NOISeq",
                "stats.SAMseq"))
  tit <- paste(com, collapse = "~")  #' plot title
  
  res <- do.call(method, list(data, cls, com, norm.method))
  #' method <-
  #'   if (is.function(method)) method
  #'   else if (is.character(method)) get(method)
  #'   else eval(method)
  #' res <- method(data, cls, com, norm.method)
  
  #' -----------------------------------------------------------------------
  #' Plot MA-plot,Volcano plot and histogram of p-values.
  tmp <- res$stats[,c("mean", "log2.fold.change","pval","padj")]
  tmp <- tmp[is.finite(tmp$"log2.fold.change"),]  #' only for panel.lmline
  
  #' MA-plot: M (log-fold change) versus A (average) plot, here showing the 
  #' features selected as differentially expressed (with a 10% false 
  #' discovery rate)
  ma.p <- 
    xyplot(log2.fold.change ~ mean, data=tmp,
           pch=16, cex=.5, 
           col=ifelse(tmp$padj < .1, "#FF000040", "black" ),
           panel = function(x, y, col,...) {
             panel.xyplot(x, y, col=col,...)     #' col="black"
             ##panel.loess(x, y, col="green",...)
             ##panel.lmline(x, y, col="red",...)   
             panel.grid(h=-1, v= 2)
             panel.abline(h = c(-2,2), col = "blue",lty =2)
           },
           xlab = "mean", ylab = "log2 fold change", 
           main = list(paste(tit, ": MA-Plot", sep=""), cex=1.0),
           scales = list(x=list(log=TRUE), y=list(log=FALSE, limits=c(-6, 6))) 
    )
  
  #' Volcano plot
  volcano.p <- 
    xyplot(-log10(pval) ~ log2.fold.change, data=tmp, 
           pch=16, cex=.5, 
           col=ifelse(abs(tmp$log2.fold.change)>2 & tmp$pval<0.05, "#FF000050", "#00000050" ),
           panel = function(x, y, col,...) {
             panel.xyplot(x, y, col=col,...)     #' col="black"
             panel.grid(h=-1, v= 2)
             panel.abline(v = c(-2,2), h=-log10(0.05),col = "blue",lty =2)
           },
           xlab = "log2(Fold Change)", ylab = "-log10(P-value)",
           main = list(paste(tit, ": Volcano Plot", sep=""), cex=1.0)
    )
  
  #' Histogram of p-values
  hist.p <- 
    histogram(~ pval, data = tmp, 
              main = list(paste(tit, ": Histogram of P-Values", sep=""), cex=1.0),
              type = "count", nint = 100)     ##"percent", "density"
  
  res$ma.p      <- ma.p
  res$volcano.p <- volcano.p
  res$hist.p    <- hist.p
  
  #' -----------------------------------------------------------------------
  #' Get reduced data according to p-values less than threshold (0.005)
  pval  <- res$stats[,"pval",drop=T]
  names(pval) <- rownames(res$stats)
  ind   <- pval <= 0.05
  
  ind   <- ifelse(is.na(ind),FALSE,ind)  #' treat NA as FALSE
  if (sum(ind) >= 1) {
    p.val <- pval[ind]
    nam   <- names(p.val)
    mat   <- res$data[,nam]
    
    dfn  <- paste(colnames(mat),"(",format(p.val, digits=3),")",sep="")
    colnames(mat) <- dfn
    
    #' MDS
    dis <- dist(mat)       #' Euclidean distance
    mds <- cmdscale(dis)
    mds <- as.data.frame(mds)
    names(mds) <- c("Coordinate_1" ,"Coordinate_2")
    mds  <- cbind(mds, cls=res$cls)
    #' MDS plot
    mds.p <- 
      xyplot(Coordinate_1 ~ Coordinate_2, data=mds,groups = cls,as.table=T, 
             xlab="Coordinate 2", ylab="Coordinate 1", 
             main=paste(tit,": MDS Plot",sep=""),
             auto.key=list(space="right"),
             par.settings = list(superpose.symbol=list(pch=rep(1:25))),
             panel=panel.elli,ep=0  
             #scales=list(cex =.75,relation="free")
      )
    
    #' PCA
    tmp  <- pca.comp(mat, scale=FALSE, pcs=1:2)
    pca  <- tmp$scores
    pca  <- cbind(pca, cls=res$cls)
    attr(pca,"varsn") <- tmp$varsn    #' attr(pcs,"varsn") #' attributes(pcs)
    #' PCA plot
    dfn <- paste(names(tmp$vars)," (", tmp$vars[names(tmp$vars)], "%)",sep = "")
    pca.p <- 
      xyplot(PC1 ~ PC2, data=pca, groups = cls, as.table=T, 
             ylab=dfn[1], xlab=dfn[2], 
             main=paste(tit,": PCA Plot",sep=""),
             auto.key=list(space="right"),
             par.settings = list(superpose.symbol=list(pch=rep(1:25))),
             panel=panel.elli,ep=0,  
             scales=list(cex =.75,relation="free")
      )
    
    #' Boxplot
    box <- cbind(mat, cls=res$cls)
    box <- melt(box, id="cls")
    box.p <- 
      bwplot(value ~ cls|variable, data = box, layout=c(3,3),
             as.table=T, pch='|',
             par.strip.text = list(cex=0.6),
             main = list(paste("Boxplot: ", tit, sep=""), cex=1.0),
             scales = list(cex=0.75,relation="free")
      )
    
    res$data.s <- mat   #' normalised data set with highly expressed variables
    res$mds    <- mds
    res$mds.p  <- mds.p
    res$pca    <- pca
    res$pca.p  <- pca.p
    res$box.p  <- box.p
  }
  return(res)
}

#' =========================================================================
#' wll-17-12-2014: Wrapper function of group stats for normalised NGS data.
#' Also the normalised data is returned.
#' Arguments: 
#'  mat    - raw count data [gene x replicate]
#'  grp    - group information
#'  nf     - normalisation factor
#'  method - data centre function
#' =========================================================================
ngs.stats <- function(mat, grp, nf, method="mean"){
  #' Normalise and transpose data
  #' nf <- nf * (1e-6*colSums(mat))       #' edgeR style  
  mat <- t(mat)/nf                
  mat <- as.data.frame(mat)
  
  #' some statistics (mean, FC, AUC, et. al)
  stats <- stats.mat(mat,grp, method=method)
  stats <- stats[,!(names(stats) %in% "pval")]
  
  return(list(stats=stats, mat=mat))
}

#' ==========================================================================
#' wll-03-12-2014: Assess feature selection of NGS: PCA plot, PLS plot and 
#'  classification.
#' ==========================================================================
ngs.cl <- function(dat.list, DF, method=c("randomForest","svm"),
                   pars=valipars(sampling = "cv",niter = 20, nreps=5, 
                                 strat=TRUE)) { 

  #' PCA and PLS plot
  pca <- pca.plot.wrap(dat.list, title=DF, par.strip.text = list(cex=0.75), 
                       scales=list(cex =.75,relation="free"),ep=0) #' free
  pls <- pls.plot.wrap(dat.list, title=DF,par.strip.text = list(cex=0.75),
                      scales=list(cex=0.75,relation="same",x=list(draw=T)))
  
  #' pdf(file = paste(PRE, tr, "pca_pls_plots.pdf", sep="_"), onefile = T)
  #' plot(pca[[1]])
  #' plot(pls[[1]])
  #' dev.off()

  #' Classification
  dn <- names(dat.list)
  #' pars <- valipars(sampling = "cv",niter = 20, nreps=5, strat=TRUE)
  #' pars <- valipars(sampling="boot",niter=1,nreps=100, strat=T)
  #' pars <- valipars(sampling="loocv")
  #' method  <- c("randomForest","svm")##,"pcalda")
  aam <- lapply(dn, function(i){
    cat("\n--data = :",i); flush.console()
    res <- aam.mcl(dat.list[[i]]$dat,dat.list[[i]]$cls,method, pars)
  })
  names(aam) <- dn
  
  #' save(aam, method, pars, file=paste(PRE,"cl.RData",sep="_"))
  #' save acc, auc and margin of classification
  ##filename  <- paste(PRE, "aam.csv", sep='_')
  ##firstline <- paste('\nACC, AUC and Margin of classification', sep='')
  ##save.tab(aam, filename, firstline)
  
  #' plot the aam results
  z <- melt(aam)
  z <- z[complete.cases(z),]  #' in case NAs
  names(z) <- c("classifier", "assessment", "value", "data")
  #' Want to remove results of pcalda?
  #' z <- subset(z, assessment!="mar")
  #' z$assessment <- factor(z$assessment)
  
  aam.p <- dotplot(factor(data, levels=rev(unique.default(data))) ~ value|assessment,
           data=z, groups = classifier, as.table=T, layout = c(length(unique(z$assessment)),1),
           par.settings = list(superpose.line = list(lty=c(1:7)),
                               superpose.symbol=list(pch=rep(1:25))),
           type="o",          #' comment this line to get original dot plot
           auto.key = list(lines=TRUE,space="bottom",columns=nlevels(z$classifier)),
           #' title="classifier"),
           xlab="", main=paste(DF," Classification with resampling", sep=":"))
  aam.p
  #' update(aam.p_fs, scale=list(x="free"),between = list(x = c(0.5, 0.5, 0), y = 0))
  #' update(aam.p_fs, scales=list(cex =.6))
  #' savePlot(filename = paste(PRE,"aam_plot", sep="_"), type = "emf")

  list(pca=pca, pls=pls, aam=aam, aam.p=aam.p)
}

#' =======================================================================
#' wll-14-08-2014: Get normalisation factors 
#' wll-04-04-2017: DESeq can be replaced by DESeq2.
#'   - method="DESeq": Low-level function to estimate size factors with
#'   robust regression.
#'   - method="TMM" is the weighted trimmed mean of M-values (to the
#'   reference) proposed by Robinson and Oshlack (2010), where the weights
#'   are from the delta method on Binomial data. If refColumn is
#'   unspecified, the library whose upper quartile is closest to the mean
#'   upper quartile is used.
#'   - method="RLE" is the scaling factor method proposed by Anders and
#'   Huber (2010). We call it "relative log expression", as median library
#'   is calculated from the geometric mean of all columns and the median
#'   ratio of each sample to the median library is taken as the scale
#'   factor.  
#'   - method="upperquartile" is the upper-quartile normalization
#'   method of Bullard et al (2010), in which the scale factors are
#'   calculated from the 75% quantile of the counts for each library, after
#'   removing genes which are zero in all libraries. This idea is
#'   generalized here to allow scaling by any quantile of the distributions.
#' =======================================================================
norm.factor <- function(data, method=c("DESeq","TMM","RLE","UQ","none"),
                        norm.data=FALSE) {
  data <- as.data.frame(data)
  method <- match.arg(method)
  nf <- 
    switch(method,
           "DESeq" = DESeq::estimateSizeFactorsForMatrix(data),        #' from DESeq
           "TMM"   = calcNormFactors(data,method="TMM"),        #' from edgeR
           "RLE"   = calcNormFactors(data,method="RLE"),        #' from edgeR
           "UQ"    = calcNormFactors(data,method="upperquartile"),#' from edgeR
           "none"  = rep(1,ncol(data))
    )
  names(nf) <- colnames(data)
  
  #' -----------------------------------------------------------------------
  #' Normalise data
  #' nf <- nf * (1e-6*colSums(mat))  #' edgeR style otherwise DESeq style
  data <- t(t(data)/nf)
  data <- as.data.frame(data)
  
  if (norm.data) {
    res <- list(nf=nf, norm.data=data)
  } else {
    res <- nf
  }
  
  return(res)  
}

#' =======================================================================
#' wll-01-12-2014: Transform count data by  two methods: 
#'  Variance stabilizing transformation (VST) and Regularized log 
#'  transformation (RLT) implemented in package DESeq2. The transformed data
#'  are also transposed for visualisation purpose.
#' =========================================================================
vst.rlt.tr <- function(data, method=c("vst", "rlt"),...){
  method <- match.arg(method)

  dat <- as.matrix(data)
  dn  <- dimnames(dat)
  dimnames(dat) <- NULL

  dat <- 
    switch(method,
           "vst" = DESeq2:::varianceStabilizingTransformation(dat,...),
           "rlt" = DESeq2:::rlog(dat,...)
           )
  dimnames(dat) <- dn
  dat <- as.data.frame(t(dat))
  return(dat)
}

#' ========================================================================
#' wll-26-11-2014:  Modify DESeq2's plotPCA: add argument select. 
#' ========================================================================
plotPCA <- function(x, intgroup="condition", ntop=500, returnData=FALSE, 
                    select, sample) {

  #' wll-26-11-2014: modification goes here.
  if (missing(select)) {
    #' calculate the variance for each gene
    rv <-  genefilter::rowVars(assay(x))    #' wll: add genefilter
    #' select the ntop genes by variance
    select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  }
    
  #' perform a PCA on the data in assay(x) for the selected genes
  if (missing(sample)){
    pca <- prcomp(t(assay(x)[select,]))
  } else {
    pca <- prcomp(t(assay(x)[select,sample]))
  }

  #' the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  #' add the intgroup factors together to create a new grouping factor
  intgroup.df <- as.data.frame(colData(x)[, intgroup, drop=FALSE])
  group <- factor(apply( intgroup.df, 1, paste, collapse=" : "))

  #' assembly the data for the plot
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, intgroup.df, names=colnames(x))

  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))
}

#' ========================================================================
#' lwc-17-07-2013: Two group stats summary.
if (F) {
  data(iris)
  sel <- c("setosa", "versicolor")
  idx <- iris[,5] %in% sel
  dat <- iris[idx,1:4]
  cls <- factor(iris[idx,5])
  #' use matrix format
  stats.mat(dat,cls, method="median")
  #' or
  mat <- melt(cbind(dat,cls))
  ddply(mat, .(variable), function(x,method){
     stats.vec(x$value, x$cls,method=method)
  }, method="median")
}

#' ========================================================================
#' lwc-26-07-2013: Group stats for matrix.
#' ========================================================================
stats.mat <- function(x,y, method="mean",test.method = "wilcox.test",...) {

  x   <- as.data.frame(x, stringsAsFactors=F)
  res <- t(sapply(x, function(i) stats.vec(i,y,method,test.method,...)))
  #' res <- t(sapply(x, function(i) stats.vec(i,y,method)))
  res <- as.data.frame(res, stringsAsFactors=FALSE)
  return(res)
}

#' ========================================================================
#' lwc-25-07-2013: Group stats for numeric vector
#' wll-11-08-2014: Add overall mean
#' ========================================================================
stats.vec <- function(x,y, method= "mean",test.method = "wilcox.test",...){

  #' overall mean
  omn <- do.call(method, list(x,na.rm=TRUE))
  names(omn) <- method
  #' group mean
  gmn <- tapply(x,y, method, na.rm=TRUE)
  names(gmn) <- paste(names(gmn), method, sep=".")

  fc        <- round(.fc(gmn[2], gmn[1]), digits=2)   
  #' lwc-23-08-2013: gmn[1] is baseline
  names(fc) <- NULL
  log2.fc   <- round(.log2.fc(fc), digits=2)
  auc       <- round(.auc(x,y), digits=2)     
  p.val     <- round(.pval(x,y,test.method=test.method,...), digits=4)
  #' p.val  <- wilcox.test(x ~ y,correct = FALSE)$"p.value"  

  res <- c(omn, gmn, fold.change=fc, log2.fold.change=log2.fc, auc=auc, pval=p.val)
  return(res)
}

#' ========================================================================
#' lwc-30-07-2013: Wrapper functions for p-values from test
#' ========================================================================
.pval <- function(x, y, test.method = "oneway.test", ...){
  test.method <- if (is.function(test.method))
      test.method
  else if (is.character(test.method))
      get(test.method)
  else eval(test.method)
  pval <- test.method(x ~ y, ...)$p.value
  return(pval)
}

#' =========================================================================
#' lwc-17-07-2013: fold change 
#' ========================================================================
.fc <- function(num, denom){
  ifelse(num >= denom, num/denom, -denom/num)
}

#' =========================================================================
#' lwc-17-07-2013: transform fold change with log2
#' =========================================================================
.log2.fc <- function(foldchange, base = 2){
  res <- ifelse(foldchange < 0, 1/-foldchange, foldchange)
  res <- log(res, base)
  res
}


#' =========================================================================
#' lwc-29-07-2013: AUC functions
#' =========================================================================
.auc <- function(stat,label,pos=levels(as.factor(label))[2]){
  label <- factor(label)
  if (nlevels(label) != 2) stop("'label' must be two categorical data")
  
  #' sort out which level is pos and convert it to "1".
  pos <- match.arg(pos, levels(label)) 
  label <- relevel(label, ref=pos)   #' put pos in the 1st position
  levels(label) <- list("0"=levels(label)[2], "1"=levels(label)[1])

  label <- label[order(stat, decreasing=T)]
  label <- as.numeric(levels(label))[as.integer(label)]

  tmp   <- cumsum(label)/sum(label)
	auc   <- mean(tmp[label==0])
  
  auc[auc < 0.5] <- 1 - auc[auc < 0.5]  
  return(auc)
}


#' =========================================================================
#' wll-17-11-2008: Get number of rejected hypotheses for several multiple
#'                 testing procedures based on Type I error rates.
#' =========================================================================
pval.reject <- function(adjp,alpha){
  adjp <- as.data.frame(adjp)
  tmp <- sapply(alpha, function(y){
    p.num <- sapply(adjp, function(x) sum(x <= y,na.rm = TRUE))
  })
  #' colnames(tmp) <- alpha
  colnames(tmp) <- paste("<=",alpha,sep="")
  tmp <- t(tmp)
  return(tmp)
}

#' ========================================================================
#' wll-26-04-2008: my version of unlist.
#' =========================================================================
un.list <- function(x, y=""){
  res <- list()
  for (i in names(x)){
    id <- if(y=="") i else paste(y,i,sep="_")
    if (is.list(x[[i]]) && !is.data.frame(x[[i]])) {
    #' Since data frame has also property of list
      tmp <- un.list(x[[i]], y=id)
      res <- c(res,tmp)
    } else {
      res[[id]] <- x[[i]]
    }
  }
  res
}

#' =========================================================================
#' wll-17-08-2014: wrapper function for shrinking list
#' Note: This function is very specific. Any recursive programming?
#'       See un.list
#' =========================================================================
sh.list <- function(x){
  tmp <- lapply(x, function(y){
    y[!sapply(y, is.null)]
  })
  tmp <- tmp[!sapply(tmp, function(y) length(y) == 0)]
}

#' =========================================================================
#' wll-29-03-2008: Compute the PCA scores and proportion of variance
#' =========================================================================
pca.comp <- function(x, scale=FALSE, pcs=1:2,...){
  pca  <- prcomp(x, scale.=scale,...)
  vars <- pca$sdev^2          #' i.e. eigenvalues/variance
  vars <- vars/sum(vars)      #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100,2)
  dfn  <- paste(names(vars),": ",vars[names(vars)],"%",sep="")
  x    <- data.frame(pca$x)
  x    <- x[,pcs]
  vars <- vars[pcs]
  dfn  <- dfn[pcs]
  return(list(scores=x,vars=vars,varsn=dfn))
}

#' =========================================================================
#' lwc-14-02-2010: Panel function for plotting ellipse  used by lattice.
#' =========================================================================
panel.elli <- function(x, y, ep=0,conf.level = 0.975, ...) {
  plot.elli <- function(x,y,...){ #' plot ellipse
    Var  <- var(cbind(x,y))
    Mean <- cbind(mean(x),mean(y))
    Elli <- ellipse(Var, centre = Mean, level = conf.level)
    #' panel.xyplot(x, y,...)
    #' panel.xyplot(Elli[,1], Elli[,2],...)
    panel.points(Elli[,1], Elli[,2],...)
  }

  panel.superpose(x,y,...)
  panel.abline(h=0, v=0,col=c("gray"), lty=2)
  if (ep == 1){
    #' overline ellipse
    plot.elli(x,y,type="l",col="red",lwd=2,...)
  } else if (ep == 2){
    #' group ellipse
    panel.superpose(x,y,..., panel.groups = plot.elli, type="l", lty=2)
  }
}

#' =========================================================================
#' wll-17-12-2014: Correct p-value with independent filtering. 
#' wll-06-01-2015: Each NGS wrapper function calls this function.
#' Note: This function is modified from some code segments from 
#'       R/Biocondunctor 'DESeq2' and 'genefilter'.
#' =========================================================================
p.adj.f <- function(pval, filter.stat, alpha=0.1, method = "BH"){
  #' Probabilities for calculating quantiles which are used as cutoff points
  #' for filtering. 
  theta  <- seq(0, 0.95, by=0.05)
  #' call genefilter's filtered_p
  padj   <- genefilter::filtered_p(filter=filter.stat, test=pval, theta=theta, 
                                   method = method)
  #' get rejection number
  rej    <- colSums(padj < alpha, na.rm = TRUE )

  #' select the maximum rejection number
  if (T) {  #' use the last tie
    ind  <- max.col(t(as.matrix(rej)), ties.method="last")
  } else {
    ind  <- which.max(rej)    #' only use the first tie
  }
  #' the final choice
  padj.f   <- padj[, ind, drop=TRUE]     

  #' get filtering cutoffs
  cutoff <- quantile(filter.stat, theta ) 
  cutoff.f <- cutoff[ind]

  list(padj.f=padj.f, cutoff.f=cutoff.f, cutoff=cutoff, 
       rej=rej, padj=cbind(raw=pval, padj))
}

#' =========================================================================
#' wll-06-01-2015: data analysis using limma + voom
#' Note: data is [gene x replicate]
#' =========================================================================
stats.voom <- function(data, cls, com, norm.method="TMM"){
  #' -----------------------------------------------------------------------
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' Construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' --------------------------------------------------------------------
  #' Prepare data for voom+limma
  design <- model.matrix(~grp)
  v   <- voom (mat, design, lib.size = colSums(mat)*nf)
  res <- lmFit(v,design)
  res <- eBayes(res)
  tab <- topTable(res,coef=ncol(design),n=nrow(mat),sort.by="none")

  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  
  
  #' ---------------------------------------------------------------------
  #' Get p-values, adjusted p-values and filtered-adjusted p-values.
  pval <- tab$P.Value 
  padj <- p.adjust(pval,method="BH")
  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  #' update results
  res$padj   <- padj
  res$padj.f <- padj.f  

  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, pval.adj)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-19-11-2014: data analysis using DESeq2
#' wll-17-12-2014: add adjusted p-values with filtering
#' wll-04-04-2017: add DESeq2:: in front of its functions.
#' Note: data is [gene x replicate]
#' =========================================================================
stats.DESeq2 <- function(data, cls, com, norm.method="TMM"){

  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' Construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' Prepare data for DESeq2
  colDat <- data.frame(condition=grp)  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                        colData = colDat,
                                        design = ~ condition)
  
  sizeFactors(dds) <- nf            #' user-defined norm factor 
  #' colData(dds)$sizeFactor <- nf  
  #' sizeFactors(dds)
  #' colData(dds)    #' here colData is a function in package GenomicRanges
  
  #' dds <- estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds,fitType = "local")  
  #' wll-04-04-2017: change fitType as "local"
  dds <- DESeq2::nbinomWaldTest(dds)     #' DESeq's default algorithm

  #' You can use the wrapper function     
  #' dds <- DESeq(dds)

  res <- results(dds)         #' it is DataFRame class in package S4Vectos
  res <- as.data.frame(res)

  #' change column name
  ind <- which(names(res) == "pvalue")
  names(res)[ind] <- "pval"

  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  
  
  #' Get p-values, adjusted p-values and filtered-adjusted p-values.
  pval   <- res$pval
  padj   <- p.adjust(pval,method="BH")
  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  #' update results
  res$padj   <- padj
  res$padj.f <- padj.f  

  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, res) #' stats <- cbind(stats, pval.adj)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-22-09-2014: Wrapper function for SAMseq
#' wll-17-12-2014: add adjusted p-values with filtering
#' Note: adjusted p-values with filtering may be not reasonable
#' =========================================================================
stats.SAMseq <- function(data, cls, com, norm.method="TMM"){

  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac         <- norm.factor(data, method=norm.method)

  #' construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]
  
  #' Call SAMseq
  set.seed(100)
  res <- SAMseq(x=mat, y=grp, resp.type= "Two class unpaired", 
                geneid=rownames(mat), genenames=rownames(mat), 
                nperms= 1000,fdr.output=1)
  #' Note: x must be count number.
  
  #' Get stats
  tmp <- res$siggenes.table
  tab <- rbind(res$siggenes.table$genes.up, res$siggenes.table$genes.lo)

  #' Get scores and q-values
  score <- rep(0, nrow(mat))     
  score[match(tab[, 1], rownames(mat))] = as.numeric(tab[, 3])
  
  fdr <- rep(1, nrow(mat))       
  fdr[match(tab[, 1], rownames(mat))] = as.numeric(tab[, 5])/100

  padj <- fdr    
  pval <- fdr    #' No p-values and therefore assign fdr to it

  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  

  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f

  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, pval.adj)

  #' order based on pval
  stats <- stats[order(stats$pval), ]

  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-15-09-2014: Wrapper function for NOISeq
#' wll-17-12-2014: add adjusted p-values with filtering
#' Note: adjusted p-values with filtering may be not reasonable
#' =========================================================================
stats.NOISeq <- function(data, cls, com, norm.method="TMM"){

  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac         <- norm.factor(data, method=norm.method)

  #' construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' NOISeq allow the normalised data being used.
  
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat   <- tmp$mat  
  
  #' Prepare data for NOISeq
  fact <- as.data.frame(grp)
  rownames(fact) <- rownames(mat)

  dat <- readData(data=t(mat), factors=fact)
  #' slotNames(dat)
  #' pData(dat)
  #' head(assayData(dat)$exprs)

  #' Call NOISeq. No normalisation and filtering. Beware that do not use
  #'  noiseq since the probability cannot be convert to FDR.
  res <- noiseqbio(dat, norm = "n", factor="grp", filter = 0, k=0.0)
  tab      <- res@results[[1]]
  tab$pval <- 1 - tab$prob
  tab$padj <- 1 - tab$prob

  #' get the adjusted p-values with filtering
  pval <- tab$pval
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  tab$padj.f <- padj.f

  #' Get p-values and FDRs
  pval.adj <- cbind(pval=tab$pval, padj=tab$padj,padj.f=tab$padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, tab)

  #' order based on pval
  stats <- stats[order(stats$pval), ]

  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-13-08-2014: Wrapper function for TSPM
#' wll-06-01-2015: minor changes
#' Note: mat is [gene x replicate]
#' =========================================================================
stats.TSPM <- function(data, cls, com, norm.method="TMM"){
  #' -----------------------------------------------------------------------
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' -----------------------------------------------------------------------
  #' call TSPM
  x0  <- rep(1, times=length(grp)) 
  res <- TSPM(mat, grp, x0, lib.size=colSums(mat)*nf) 
  #' user-defined norm factor
  
  #' names(res)  
  #' [1] "log.fold.change"     "pvalues"  "index.over.disp"    
  #' [4] "index.not.over.disp" "padj" 
  #' wll-24-04-14: Since we cannot convert res into data frame,only 
  #'  "pvalues" and "padj" are returned.

  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  

  #' ---------------------------------------------------------------------
  #' Get p-values and FDRs
  pval <- res$pvalues
  padj <- res$padj
  #' padj   <- p.adjust(pval,method="BH")
  #' wll-06-01-2015: two versions of padj are different.
  
  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f

  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)

  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  #' stats <- cbind(stats, res)
  stats <- cbind(stats, pval.adj)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-13-08-2014: Wrapper function for edgeR
#' wll-06-01-2015: minor changes
#' Note: mat is [gene x replicate]
#' To-Do: 
#'  1.) update 'stats.edgeR' with 'glmFit' and 'glmLRT'. 
#'  2.) consider the design matrices and contrasts for multi-factor test
#' =========================================================================
stats.edgeR <- function(data, cls, com, norm.method="TMM"){

  #' -----------------------------------------------------------------------
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' ------------------------------------------------------------------------
  #' Prepare data for edgeR
  dge <- DGEList(counts=mat, group=grp) 
  dge$samples$norm.factors <- nf      #' user-defined normalisation factor
  #' dge <- calcNormFactors(dge)

  dge <- estimateCommonDisp(dge)
  dge <- estimateTagwiseDisp(dge)
  
  #' ------------------------------------------------------------------------
  #' wl-27-07-2019, Sat: check 'glmLRT'
  res <- exactTest(dge)   #' default: pair=1:2
  tab <- res$table
  names(tab)[3] <- "pval"
  
  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  
  
  #' ---------------------------------------------------------------------
  #' Get p-values, adjusted p-values and filtered-adjusted p-values.
  pval <- tab$pval
  padj <- p.adjust(pval, method="BH")
  #' wll-04-08-14: stats (ordered by pval) can be obtained directly by
  #'  topTags(res, n=nrow(res$table))

  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  #' update results
  tab$padj   <- padj
  tab$padj.f <- padj.f  

  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, tab)
  #' stats <- cbind(stats, pval.adj)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  #' wl-27-07-2019, Sat: return `edgeR` dge?
  return(res)
}

#' =========================================================================
#' wll-13-08-2014: data analysis using DESeq
#' wll-06-01-2015: minor changes
#' wll-04-04-2017: add DESeq:: to avoid conflict with DESeq2 
#' Note: mat is [gene x replicate]
#' =========================================================================
stats.DESeq <- function(data, cls, com, norm.method="TMM"){
  #' -----------------------------------------------------------------------
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' Construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' --------------------------------------------------------------------
  #' Prepare data for DESeq
  cds <- DESeq::newCountDataSet(mat, grp) 
  sizeFactors(cds) <- nf            #' user-defined norm factor
  #' cds <- DESeq::estimateSizeFactors(cds)
  cds <- DESeq::estimateDispersions(cds,fitType = c("local"))   
  #' default method is "pooled"

  #' Call DESeq
  res <- DESeq::nbinomTest(cds, condA=levels(grp)[1], condB=levels(grp)[2])
  
  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  

  #' ---------------------------------------------------------------------
  #' Get p-values, adjusted p-values and filtered-adjusted p-values.
  pval <- res$pval
  padj <- res$padj
  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  #' update results
  res$padj.f <- padj.f  

  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, res)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-02-09-2014: data analysis using NBPSeq
#' wll-06-01-2015: minor changes
#' Note: mat is [gene x replicate]
#' =========================================================================
stats.NBPSeq <- function(data, cls, com, norm.method="TMM"){
  #' -----------------------------------------------------------------------
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' Construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' --------------------------------------------------------------------
  #' Call NBPSeq
  res <- nbp.test(as.matrix(mat), grp, grp1=levels(grp)[1], grp2=levels(grp)[2], 
                  norm.factors = nf)
  #' Show boxplots, MA-plot, mean-variance plot and mean-dispersion plot
  #' par(mfrow=c(3,2));
  #' plot(res);
  
  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  

  #' ---------------------------------------------------------------------
  #' Get p-values, adjusted p-values and filtered-adjusted p-values.
  pval <- res$p.values
  #' padj <- res$q.values    #' wll-15-09-2014: qvalues is not equal to FDR.
  padj <- p.adjust(pval, method="BH")
  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, pval.adj)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-03-09-2014: data analysis using EBSeq
#'   Original author's claim: In EBSeq 1.3.3, the default setting of EBTest
#'   function will remove low expressed genes (genes whose 75th quantile of
#'   normalized counts is less than 10) before identifying DE genes. These
#'   two thresholds can be changed in EBTest function. We found that low
#'   expressed genes are more easily to be affected by noises. Removing
#'   these genes prior to downstream analyses can improve the model fitting
#'   and reduce impacts of noisy genes (e.g. genes with outliers).
#' wll comment: Totally agree. But any filtering should be done before
#'   calling the differential analysis function
#' wll-06-01-2015: minor changes. padj.f may not be reasonable.
#' =========================================================================
stats.EBSeq <- function(data, cls, com, norm.method="TMM"){
  #' -----------------------------------------------------------------------
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac    <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' Construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' --------------------------------------------------------------------
  #' Call EBSeq
  res <- EBTest(Data=as.matrix(mat), Conditions=grp, sizeFactors=nf, maxround=5,
                QtrmCut=-1)  
  #' Note: 1.) Data must be matrix
  #'       2.) Do NOT filter genes by quantile, so QtrmCut=-1. The filtering 
  #'           should be done outside this function.
  pp  <- GetPPMat(res)
  ppde <- 1 - pp[, "PPDE"]   #' Is this FDR?

  #' Get FDR. Modified from R scripts of A comparison of methods for differential
  #' expression analysis of RNA-seq data by Charlotte Soneson and Mauro Delorenzi
  #' BMC Bioinformatics 2013
  fdr <- sapply(ppde,function(x){
     mean(ppde[which(ppde <= x)])
  })
  
  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat <- tmp$mat  

  #' ---------------------------------------------------------------------
  #' Get p-values, adjusted p-values and filtered-adjusted p-values.
  pval   <- fdr
  padj   <- fdr
  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)      #' will pass names of p-values
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f
  
  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, pval.adj)
  
  #' order based on pval
  stats <- stats[order(stats$pval), ]
  
  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-13-08-2014: Wrapper function for Wilcoxon test
#' wll-06-01-2015: minor changes. 
#' Note: mat is [gene x replicate]
#' =========================================================================
stats.Wilcox <- function(data, cls, com, norm.method="TMM"){
  #' Get normalisation factors
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  fac         <- norm.factor(data, method=norm.method)

  #' -----------------------------------------------------------------------
  #' construct binary comparison info
  tmp <- paste(com, collapse = "|")
  ind <- grepl(tmp, cls)
  mat <- data[,ind,drop=FALSE] 
  nf  <- fac[ind]  
  grp <- cls[ind, drop=TRUE]

  #' -----------------------------------------------------------------------
  #' Wilcox allow the normalised data being used.
  #' --------------------------------------------------------------------
  #' get group stats and transposed-normalised data
  tmp   <- ngs.stats(mat, grp, nf, method="mean")
  stats <- tmp$stats
  #' transposed and normalised data
  mat   <- tmp$mat  
  
  #' -----------------------------------------------------------------------
  #' call Wilconxon Test
  pval <- sapply(as.data.frame(mat), function(x){
    .pval(x, grp, test.method = "wilcox.test", correct = FALSE)
  })
  #' wll-01-07-2014: 1. correct should be FLASE;
  #'                 2. p-values can be extracted from stats.
  #' adjust p-values (BH: Benjamini¨CHochberg. i.e. fdr)
  padj <- p.adjust(pval, method="BH")

  #' get the adjusted p-values with filtering
  names(pval) <- colnames(mat)
  padjf  <- p.adj.f(pval, filter.stat=stats$mean, alpha=0.1, method = "BH")
  padj.f <- padjf$padj.f

  #' ------------------------------------------------------------------
  pval.adj <- cbind(pval=pval, padj=padj, padj.f=padj.f)
  pval.adj <- data.frame(pval.adj)  
  rownames(pval.adj) <- colnames(mat)
  
  #' get reject number
  tmp     <- pval.adj[order(pval.adj[,1],na.last=T),]  
  rej.num <- pval.reject(tmp, alpha=c(0.01, 0.05, 0.1)) 

  #' merge adjusted p-values
  stats <- cbind(stats, pval.adj)

  #' order based on pval
  stats <- stats[order(stats$pval), ]

  res <- list(stats=stats, rej.num=rej.num, data=mat, cls=grp, padjf=padjf)
  return(res)
}

#' =========================================================================
#' wll-23-09-2014: Replace the original samr.estimate.depth. Beware that 
#'  argument must be norm.method for keeping consistent with the original one
#'  and passing the argument outside. This function should also been assigned 
#'  in the name space.
#' =========================================================================
samr.estimate.depth <- function(x,method=norm.method) {
  norm.method <- match.arg(norm.method, c("DESeq","TMM","RLE","UQ","none"))
  print(norm.method)
  fac         <- norm.factor(x, method=norm.method)
  return(fac)
}
#' override/replace a function in a package namespace
assignInNamespace("samr.estimate.depth",samr.estimate.depth, ns="samr")

#' =========================================================================
#' wll-11-12-2014: Wrapper function to get OTUs table
#' wll-03-02-2015: Fixed a bug for NAs
#' to-do: frequency of each OTU identified.
#' =========================================================================
otu.wrapper <- function(res, pval="padj", alpha=0.1){
  #' get otus
  otu <- lapply(res, function(x) sig.var(x$res,pval, alpha))
  otu <- otu[!is.na(otu)]  
  #' wll-03-02-2015:  remove NAs before melting. should be a efficient and 
  #' general method to shrink list. My version of shrink.list in mt has 
  #' problem since recent version of R. 
  
  #' get contingency tables
  tmp <- melt(otu)  
  #' remove NAs if any
  #' tmp <- tmp[complete.cases(tmp),]
  
  otu_num <- cast(tmp,   L1 ~ ., value="value")
  names(otu_num) <- c("design", "otu_number")
  
  otu_tab <- cast(tmp,   value ~ L1)
  names(otu_tab)[1] <- "otu"
  
  #' summary of OTUs table
  otu_all <- do.call("c", otu)
  otu_tab_summ <- table(otu_all)
  otu_freq <- as.data.frame(table(otu_tab_summ))
  
  list(otu_tab=otu_tab,otu_num=otu_num, otu_freq=otu_freq)
}

#' =========================================================================
#' wll-10-12-2014: get significant features (i.e. genes or OTUs)
#' =========================================================================
sig.var <- function(x, pval="padj", alpha=0.1) {
  #' order stats
  stats <- x[order(x[[pval]]),]
  idx   <- which(stats[[pval]] <= alpha)
  if (length(idx) != 0) {
    tmp   <- stats[idx, pval, drop=F]
    var   <- rownames(tmp)
    return(var)
  } else {
    return(NA)
  }
}

#' =========================================================================
#' wll-09-12-2014: Summary function for NGS vector data
#' =========================================================================
vec.summ.ngs <- function(x){  
  res <- c(N = sum(!is.na(x)), 
           N.nonzero = sum(x > 0),
           Sum = sum(x,na.rm=T),
           Min = min(x,na.rm=T), 
           Mean = mean(x, na.rm=T), 
           Median = median(x, na.rm=T),
           Max = max(x, na.rm=T), 
           Std = sd(x, na.rm=T))
  res <- round(res, digits=3)
  return(res)
}

#' =========================================================================
#' wll-02-07-2014: TSPM algorithm
#' 	R code for the paper by Paul L. Auer and R.W. Doerge:
#' 	"A Two-Stage Poisson Model for Testing RNA-Seq Data"
#' 	Example:
#' Input:
#'  counts:   a matrix of RNA-Seq gene counts (genes are rows, samples are 
#'            columns)
#'  x1: 		  a vector of treatment group factors (under the alternative 
#'            hypothesis)
#'  x0: 		  a vector of treatment group factors (under the null hypothesis)
#'  lib.size: a vector of RNA-Seq library sizes. This could simply be obtained
#'          	by specifying lib.size <- apply(counts,2,sum). It may also be 
#'            any other	appropriate scaling factor.
#'  alpha.wh:	the significance threshold to use for deciding whether a gene 
#'            is overdispersed. Defaults to 0.05.
#' Output:
#'  log.fold.change:		 a vector containing the estimated log fold changes 
#'                       for each gene
#'  pvalues: 			       a vector containing the raw p-values testing 
#'                       differential expression for each gene.
#'  index.over.disp: 	   a vector of integer values containing the indices of 
#'                       the over-dispersed genes.
#'  index.not.over.disp: a vector of integer values containing the indices of 
#'                       the non-over-dispersed genes.
#'  padj:			           a vector containing the p-values after adjusting for 
#'                       multiple testing using the method of 
#'                       Benjamini-Hochberg
#' ========================================================================
TSPM <- function(counts, x1, x0, lib.size, alpha.wh=0.05){

  #' Initializing model parameters
  n <- dim(counts)[1]
  per.gene.disp <- NULL
  LRT <- NULL
  score.test <- NULL
  LFC <- NULL
  
  #' Fitting the GLMs for each gene
  for(i in 1:n){
    #' Fit full and reduced models
    model.1 <- 
      glm(as.numeric(counts[i,]) ~ x1, offset=log(lib.size), family=poisson)
    model.0 <- 
      glm(as.numeric(counts[i,]) ~ x0, offset=log(lib.size), family=poisson)
    
    #' Obtain diagonals of Hat matrix from the full model fit
    hats <- hatvalues(model.1)
    
    #' Obtain Pearson overdispersion estimate
    per.gene.disp[i] <- 
      sum(residuals(model.1, type="pearson")^2)/model.1$df.residual
    
    #' Obtain Likelihood ratio statistic
    LRT[i] <- deviance(model.0) - deviance(model.1)
    
    #' Obtain score test statistic
    score.test[i] <- 
      1/(2*length(counts[i,])) * sum(residuals(model.1, type="pearson")^2 - ((counts[i,] - hats*model.1$fitted.values)/model.1$fitted.values))^2
    
    #' Obtain the estimated log fold change
    #' lwc-03-07-2014: Need to check it out.
    LFC[i] <- -model.1$coef[2]
  }
  
  #' Initialize parameters for Working-Hotelling bands around the score TSs 
  qchi <- qchisq(df=1, (1:n-0.5)/n)
  MSE <- 2
  UL <- NULL
  
  #' Obtain the upper boundary of the WH bands 
  xbar   <- mean(qchi)
  bottom <- sum((qchi-xbar)^2)
  top    <- (qchi-xbar)^2
  s      <- sqrt(MSE*(1/n) + (top/bottom))
  W      <- sqrt(2*qf(df1=1, df2=n-1, p=1-(alpha.wh/n)))
  UL     <- pmax(qchi + W*s,1)
  
  #' Obtain the indices of the over-dispersed and not-over-dispersed genes, 
  #' respectively 
  cutoff <- min(which(sort(score.test)-UL > 0))
  temp <- cutoff-1 + seq(cutoff:length(score.test))
  over.disp <- which(score.test %in% sort(score.test)[temp])
  not.over.disp <- setdiff(1:length(score.test), over.disp)
  
  #' Compute p-values
  p.f   <- pf(LRT[over.disp]/per.gene.disp[over.disp], df1=1, #' lwc: The F Distribution
              df2=model.1$df.residual, lower.tail=FALSE)
  p.chi <- pchisq(LRT[not.over.disp], df=1, lower.tail=FALSE)
  p     <- NULL
  p[over.disp] <- p.f
  p[not.over.disp] <- p.chi
  
  #' Adjust the p-values using the B-H method
  p.bh.f <- p.adjust(p.f, method="BH")
  p.bh.chi <- p.adjust(p.chi, method="BH")
  final.p.bh.tagwise <- NULL
  final.p.bh.tagwise[over.disp] <- p.bh.f
  final.p.bh.tagwise[not.over.disp] <- p.bh.chi
  
  #' Output
  list(log.fold.change=LFC, pvalues=p, index.over.disp=over.disp, 
       index.not.over.disp=not.over.disp, padj=final.p.bh.tagwise)
}


#' ==========================================================================
#' Example of TSPM
#' ==========================================================================
if (F) {
  counts <- matrix(0, nrow=1000, ncol=10)
  for(i in 1:1000){
    lambda <- rpois(n=1, lambda=10)
    counts[i,] <- rpois(n=10, lambda=lambda)
  }
  
  x1 <- gl(n=2, k=5, labels=c("T", "C"))  #' [1] T T T T T C C C C C
  x0 <- rep(1, times=10)                  #' [1] 1 1 1 1 1 1 1 1 1 1 (lwc:?)
  lib.size <- apply(counts,2,sum)
  
  res <- TSPM(counts, x1, x0, lib.size)
  
  #' extreive results
  pval <- res$pvalues
  padj <- res$padj
  
  tmp <- cbind(pval=pval, padj=padj)
  tmp <- tmp[order(tmp[,1]),]
}

#' ==========================================================================
#' lwc-06-06-2013:  Taken from Multiple graphs on one page (ggplot2), Cookbook 
#' for R. http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' Multiple plot function
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot 
#' objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' ==========================================================================
multiplot <- function(p.list, plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  #' Make a list from the ... arguments and plotlist
  #' plots <- c(list(...), plotlist)
  plots <- c(p.list, plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    #' Make the panel
    #' ncol: Number of columns of plots
    #' nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    #' Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    #' Make each plot, in the correct location
    for (i in 1:numPlots) {
      #' Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' =======================================================================
#' Function List
#' =======================================================================
#' (01). otu.wrapper 
#' (02). sig.var 
#' (03). vec.summ.ngs 
#' (04). ngs.cl 
#' (05). vst.rlt.tr 
#' (06). plotPCA 
#' (07). stats.mat 
#' (08). stats.vec 
#' (09). .pval
#' (10). .fc
#' (11). .log2.fc
#' (12). .auc
#' (13). pval.reject 
#' (14). un.list 
#' (15). sh.list 
#' (16). pca.comp 
#' (17). panel.elli 
#' (19). multiplot 
#' (23). ngs.stats 
#' (24). p.adj.f 
#' (25). stats.voom 
#' (26). stats.DESeq2 
#' (27). stats.SAMseq 
#' (28). stats.NOISeq 
#' (29). stats.TSPM 
#' (30). stats.edgeR 
#' (31). stats.DESeq 
#' (32). stats.NBPSeq 
#' (33). stats.EBSeq 
#' (34). stats.Wilcox 
#' (35). TSPM 
#' (36). samr.estimate.depth 
#' (37). ngs 
#' (38). norm.factor 
