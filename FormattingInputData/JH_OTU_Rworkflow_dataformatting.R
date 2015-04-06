env=read.table("environ.txt", header=TRUE, sep="\t", row.names=1)
library(vegan)

#standardize env2 (z-score) and transpose matrix
env2=decostand(env,method="standardize", MARGIN=2)
env2=t(env2)

#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("JH_OTUs.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "JH_OTUs_nosigs_env.txt", sep="\t", quote=FALSE)

