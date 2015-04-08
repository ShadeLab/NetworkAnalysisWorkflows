#Determine P values based on alpha= 0.6 tables from MINE
mine_fp="JH_OTUs_nosigs2.csv,allpairs,cv=0.0,B=n^0.6,Results.csv"
pvalue_fp="n=25,alpha=0.6_forR.csv"

#download the pvalue table, delete the first 5-6 rows and also the 5-6 "invisible" last rows. colnames = "MIC", "P", "conf"

MINE_p.f=function(mine_fp,pvalue_fp){
    mine=read.csv(mine_fp, header=TRUE, check.names=FALSE)
    pv=read.csv(pvalue_fp, header=TRUE,check.names=FALSE)
    p=NULL
    for (i in 1:nrow(mine)){
        mic=mine[i, "MIC (strength)"]
        pvMIC=pv[,"MIC"]
        estp=mic > pvMIC
        estp.0=pv[estp,]
        estp.1=estp.0[1,"P"]
        p[i]=estp.1
    }
    
    out=cbind(mine,p)
    colnames(out)=c(colnames(mine), "p_value_uncorrected")
    write.table(p, "MINE_pvalue.txt", quote=FALSE, sep="\t", row.names=FALSE)
    return(p)
    return(out)
}

#run the funciton
MINE_p.f(mine_fp, pvalue_fp)

#Correct p values using the false discovery rate
#download and install multtest package from bioconductor
#add link here
fdr.f=function(pvalue_fp){
    library(multtest)
    pv=read.csv(pvalue_fp, header=TRUE)
    p=pv[,1]
    out=mt.rawp2adjp(p,proc=c("BH", "BY"), alpha=0.05, na.rm=TRUE)
    #return(out)
    bh=out$adjp[(out$adjp[,"BH"]>0.054),]
    by=out$adjp[(out$adjp[,"BY"]>0.054),]
    print(grep(bh[1,"BH"],out$adjp[,"BH"]))
    print(grep(by[1,"BY"],out$adjp[,"BY"]))
    write.table(out$adjp, "AdjP.txt",quote=FALSE, sep="\t", row.names=FALSE)
}

pvalue_fp="MINE_pvalue.txt"
fdr.f(pvalue_fp)

#Append adjusted pvalue to original file
mine_fp="JH_OTUs_nosigs2.csv,allpairs,cv=0.0,B=n^0.6,Results.csv"
adjp_fp="AdjP.txt"

mine=read.csv(mine_fp, header=TRUE, check.names=FALSE)
adjp=read.table(adjp_fp, header=TRUE, check.names=FALSE)
combined=cbind(mine,adjp)
write.table(combined, "JH_OTUs_nosigs2.csv,allpairs,cv=0.0,B=n^0.6,Results,adjp.txt", quote=FALSE, sep="\t")

#Convert MINE CSV output to Cytoscape format
mine2Cytoscape.f=function(mine_fp, name){
    data=read.table(mine_fp, header=TRUE, check.names=FALSE, sep="\t")
    data2=data[(2*is.na(data[,"BH"])==0),]
    write.table(data2,paste(name,"_mine2Cyto.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
}
mine_fp="JH_OTUs_nosigs2.csv,allpairs,cv=0.0,B=n^0.6,Results,adjp.txt"
mine2Cytoscape.f(mine_fp, name="JH_OTUs_nosigs2")
