env=read.table("Cen_Map_3.txt", header=TRUE, sep="\t", row.names=1)
library(vegan)

#standardize env2 (z-score) and transpose matrix
env2=decostand(env,method="standardize", MARGIN=2)
env2=t(env2)

# For Geobacillus
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("Geobacillus_only.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Geobacillus_only_no_sigs.txt", sep="\t", quote=FALSE)


# For Crenarchaeota
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("Crenarchaeota_only.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Crenarchaeota_only_no_sigs.txt", sep="\t", quote=FALSE)



# For L2
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("otu_table_mc2_w_tax_even455113_L2.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Cen_OTUs_nosigs_env_L2.txt", sep="\t", quote=FALSE)

# For L3
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("otu_table_mc2_w_tax_even455113_L3.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Cen_OTUs_nosigs_env_L3.txt", sep="\t", quote=FALSE)

# For L4
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("otu_table_mc2_w_tax_even455113_L4.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Cen_OTUs_nosigs_env_L4.txt", sep="\t", quote=FALSE)

# For L5
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("otu_table_mc2_w_tax_even455113_L5.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Cen_OTUs_nosigs_env_L5.txt", sep="\t", quote=FALSE)

# For L6
#reduce the dataset by removing singleton OTUs and complete absences
data=read.table("otu_table_mc2_w_tax_even455113_L6.txt", header=TRUE, sep="\t", row.names=1)
c=colSums(data)
data.pa=1*data>0
r=rowSums(data.pa)

data.nosigs=data[r>1,]

#append the environ data to the dataframe of otus (no sigs)
combined=rbind(data.nosigs,env2)

write.table(combined, "Cen_OTUs_nosigs_env_L6.txt", sep="\t", quote=FALSE)

