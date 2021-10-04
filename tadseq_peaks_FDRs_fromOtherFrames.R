#### For pvl.txt files create p-val fdr using all the other strands as the negative controls
#### this was put into an rscript that is framesFDR.R after tested to run in parallel for all desired samples

args <- commandArgs()

exp <- paste(sub('--Sample=', '', args[grep('--Sample=', args)]))
Cores <- sub('--Cores=', '', args[grep('--Cores=', args)])
Cores <- as.numeric(Cores)

## 1 core is the R process
Cores <- Cores-1

library(parallel)

#exps <- gsub("_0\\+.peaks.bed$", "", list.files("peaks/",pattern="0\\+.peaks.bed$",full.names=F))

Dir <- "./peaks/fdrPvals/"
if(!(file.exists(Dir))) {
    dir.create(Dir,FALSE,TRUE)  
}

fdrs <- function(exp){
    posBed <- read.table(paste("peaks/",exp,"_0+.peaks.bed",sep=""))
    for(x in 0:2){ 
        if(file.info(paste("peaks/",exp,"_",x,"+.peaks.txt",sep=""))$size != 0){
            PLUS <- read.table(paste("peaks/",exp,"_",x,"+.peaks.txt",sep=""))
            assign(paste0("plus",x),PLUS)
        }else{
            PLUS <- data.frame(V1=NA, V2=NA, V3=NA, V4=NA, V6=NA, V7=NA, V8=NA, V9=NA)
            assign(paste0("plus",x),PLUS)
        }   
        if(file.info(paste("peaks/",exp,"_",x,"-.peaks.txt",sep=""))$size != 0){
            MINUS <- read.table(paste("peaks/",exp,"_",x,"-.peaks.txt",sep=""))
            assign(paste0("minus",x),MINUS)
        }else{
            MINUS <- data.frame(V1=NA, V2=NA, V3=NA, V4=NA, V6=NA, V7=NA, V8=NA, V9=NA)
            assign(paste0("minus",x),MINUS)
        }
    }
    pos <- plus0[,"V9"]
    neg <- c(plus1[,"V9"]
            ,plus2[,"V9"]
            ,minus0[,"V9"]
            ,minus1[,"V9"]
            ,minus2[,"V9"]
             )
    neg <- neg[complete.cases(neg)]
    ##make a df to populate for the uniq pvals
    sub <- data.frame(pos=sort(unique(pos), decreasing=F))
    sub$fdr <- 0
    ## calc fdrs and output into a dataframe
    fdr <- as.data.frame(do.call(rbind,mclapply(pos, function(pval){ 
        cat("pval:", pval, sep="\n")
        ## select df with a pval >= pval
        posSum <- sum(pos <= pval)
        negSum <- sum(neg <= pval)/5
        fdr <- 0
        if (negSum>0){
            fdr=round(negSum/posSum*100,2)
        }     
        return(c(pval,posSum,negSum, fdr))
    }, mc.cores=Cores, mc.preschedule=TRUE, mc.silent=TRUE)
    ))
    names(fdr) <- c("pval", "posSum", "negSum", "fdr")
    stopifnot(fdr$pval==plus0$V9)
    ## combine df with fdrSum and fdrMax
    all <- cbind(plus0,fdr[,2:4])
    iv <- match(plus0$V2, posBed$V8)
    df <- cbind(posBed[iv,1:4],all[,2:ncol(all)])
    stopifnot(fdr$pval==df$V9)
    write.table(df, file=paste(Dir, exp, "_0+_peaksfdr",".txt", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
}

invisible(fdrs(exp))
#invisible(lapply(as.list(exps),fdrs))
