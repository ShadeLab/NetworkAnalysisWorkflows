#getCRANmirrors(all = FALSE, local.only = FALSE)
options(repos=structure(c(CRAN="http://ftp.ussg.iu.edu/CRAN/")))

#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")

#biocLite("bioDist")
#biocLite("impute")
#install.packages("WGCNA")

#biocLite("preprocessCore")
#biocLite("Rsamtools")

update.packages(checkBuilt=TRUE, ask=FALSE)

library("microbenchmark")
library("snow")

#library("Matrix")
library("bbmle")
library("MASS")
library("pscl")

#library("ape")
#library("seqRFLP")
#library("Rsamtools")
#library("phylotools")
#library("boot")
#library("phyloseq")

library("flashClust")
library("WGCNA")
#library("bioDist")

library("ggplot2")
library("scales")
library("grid")
#library("ggmap")

library("vegan")
library("reshape2")
library("plyr")
#library("abind")
library("stringr")
library("gridExtra")
library("igraph")

#packageVersion("phyloseq")


### GRAPH NETWORK WITH IGRAPH ###
# load OTU abundance table. make sure value are raw counts, NOT relative abundances
otu_table <_ #XXX

# optional: subset your otu_table to exclude rare OTUs

# get the names of your OTUs into sub.otus variable
sub.otus <- as.character(rownames(otu_table))

# read in correlation matrix
setwd("~/Dropbox/work/bioinformatics/sparcc/results/summer.3.1618.atleast80")

corr = read.table("data.sparcc.txt", header=TRUE, sep="\t")
corr.mat <- subset(corr, OTU_id %in% sub.otus, select=sub.otus)
	# subset correlation matrix to only include OTUs from the otu_table
	# OTU_id is the value in the upper left cell i.e. the name of leftmost column

# read in p-values matrix
pvals = read.table("pvals_two_sided.txt", header=TRUE, sep="\t")
pvals.mat <- subset(pvals, OTU_id %in% sub.otus, select=sub.otus)

# define p-value threshold and assign zero to corr values below threshold
threshold <- (pvals.mat < 0.0001) * 1
sig.mat <- as.matrix(corr.mat * threshold)

# optional: define correlation threshold and assign zeros to corrs below threshold
#sig.mat <- sig.mat*(sig.mat>=0.5)

# optional: positive only matrix
sig.mat[sig.mat < 0] <- 0

# create graph using igraph
sparcc.graph <- graph.adjacency(sig.mat, weighted=TRUE, mode="undirected", diag=FALSE)

V(sparcc.graph)
E(sparcc.graph)

# optional: change names of verteces
names <- V(sparcc.graph)$name
sample.split <- strsplit(names,'Otu',fixed=TRUE)
sample.split.temp <- ldply(sample.split)
sample.split.temp$V2 <- sprintf("%03d", as.numeric(sample.split.temp$V2))
names.num <- sample.split.temp$V2
names <- paste("Otu", names.num, sep="")
V(sparcc.graph)$name <- names

# add color to edges
# for signed matrix
#color.mat <- sig.mat
#color.mat[color.mat > 0] <- "lightcoral"
#color.mat[color.mat < 0 & color.mat != "lightcoral"] = "lightslateblue"
# vectorize color matrix and take upper triangle
#edge.color.sym <- as.vector(t(color.mat)[as.vector(t(color.mat))!=0])
#upper.in <- length(edge.color.sym)/2
#edge.color <- edge.color.sym[1:upper.in]
# add color to edges
#E(sparcc.graph)$color <- edge.color

# for positive colors grayscale ramp
# transform data into 0:1 range with min=0 & max=1
range01 <- function(x)(x-min(x))/diff(range(x))

cRamp <- function(x){
  cols <- colorRamp(gray.colors(length(V(sparcc.graph))))(range01(exp(range01(x))))
  apply(cols, 1, function(xt)rgb(xt[1], xt[2], xt[3], maxColorValue=255))
}

E(sparcc.graph)$color <- cRamp(E(sparcc.graph)$weight)

# add abundance (weight) values to vertices (OTUs)
V(sparcc.graph)$weight <- taxa_sums(y.sum)

	# optional: display the graph
	# plot.igraph(sparcc.graph, layout=layout.fruchterman.reingold,
		# edge.color=E(sparcc.graph)$color,
		# vertex.size=exp(range01(V(sparcc.graph)$weight)),
		# vertex.color="black",
		# vertex.label=NA)

# remove unconnected vertxes: subset sparcc.graph
vert.zero <- V(sparcc.graph)[degree(sparcc.graph) < 1] # identify unconnected vertexes
sparcc.graph.sub <- delete.vertices(sparcc.graph, vert.zero) # remove unconnected vertexes

## Community Calculation ##
# sparcc.comm.between.unweighted <- edge.betweenness.community(sparcc.graph.sub, weights=NULL)

sparcc.comm.between.weighted <- edge.betweenness.community(sparcc.graph.sub, weights=E(sparcc.graph.sub)$weight)
	# weights must be positive
	# weighted edges: TRUE
	# weighted nodes: FALSE
	# pro: accurate, considered golden standard
	# con: slow

# sparcc.comm.fast <- fastgreedy.community(sparcc.graph.sub)
	# # weights must be positive
	# # weighted edges: TRUE
	# # weighted nodes: FALSE
	# # pro: fast
	# # con: poor resolution i.e. struggles to identify small communities

# sparcc.comm.walk <- walktrap.community(sparcc.graph.sub)
	# # notes: struggles with small communities
	# # weighted edges: TRUE
	# # weighted nodes: FALSE
	# # pro: fast, slightly better than fastgreedy; uses random walker: makes few assumptions about network topography
	# # con: slightly slower than fastgreedy

# sparcc.comm.prop <- label.propagation.community(sparcc.graph.sub)
	 # # weights must be positive
	 # # weighted edges: TRUE
	 # # weighted nodes: FALSE
	 # # pro: very fast
	 # # con: based on randomization procedure and therefore has to be bootstrapped, which eliminates the speed the advantage

# sparcc.comm.multi <- multilevel.community(sparcc.graph.sub)
	# # weights must be positive

# sparcc.comm.lead <- leading.eigenvector.community(sparcc.graph.sub)
	# # weights must be positive
	# # weighted edges: FALSE
	# # weighted nodes: FALSE
	# # pro: fast with results close to betweenness "gold standard"
	# # con: no weighted edges

# sparcc.comm.info.e.na.v.na <- infomap.community(sparcc.graph.sub,
	# e.weights=NA,
	# v.weights=NA)
# sparcc.comm.info.e.weight.v.na <- infomap.community(sparcc.graph.sub,
	# e.weights=E(sparcc.graph.sub)$weight,
	# v.weights=NA)
# sparcc.comm.info.e.weight.v.weight <- infomap.community(sparcc.graph.sub,
	# e.weights=E(sparcc.graph.sub)$weight,
	# v.weights=V(sparcc.graph.sub)$weight)
	# # weights must be positive
	# # weighted edges: FALSE
	# # weighted nodes: FALSE

# # slow (do not use if graph has over a couple of hundred vertices)
# sparcc.comm.opt <- optimal.community(sparcc.graph.sub)
	# # need connected graph
	# #sparcc.comm.spin <- spinglass.community(sparcc.graph.sub)

# sizes(sparcc.comm.between.weighted)
# sizes(sparcc.comm.fast)
# sizes(sparcc.comm.walk)
# sizes(sparcc.comm.prop)
# sizes(sparcc.comm.multi)
# sizes(sparcc.comm.lead)
# sizes(sparcc.comm.info.unweighted)
# sizes(sparcc.comm.info.weighted)

# # visualize communities
# plots <- list(
	# sparcc.comm.between.unweighted,
	# sparcc.comm.between.weighted,
	# sparcc.comm.info.e.na.v.na,
	# sparcc.comm.info.e.weight.v.na,
	# sparcc.comm.info.e.weight.v.weight
	# )

# titles <- list(
	# "sparcc.comm.between.unweighted",
	# "sparcc.comm.between.weighted",
	# "sparcc.comm.info.e.na.v.na",
	# "sparcc.comm.info.e.weight.v.na",
	# "sparcc.comm.info.e.weight.v.weight"
	# )

 #l.f <- layout.fruchterman.reingold(sparcc.graph.sub)
 #l.a <- layout.auto(sparcc.graph.sub)

# par(mfrow=c(2,3))

# plot.igraph(sparcc.graph.sub,
	# main="sparcc.graph",
	# layout=l.a,
	# edge.color=E(sparcc.graph.sub)$color,
	# vertex.size=exp(range01(V(sparcc.graph.sub)$weight))*3,
	# vertex.color="black",
	# vertex.label=NA
# )

# for(i in 1:5){	
	# plot(plots[[i]], sparcc.graph.sub,
		# main=titles[[i]],
		# layout=l.a,
		# vertex.size=exp(range01(V(sparcc.graph.sub)$weight))*3,
		# vertex.label=NA,
		# edge.color=E(sparcc.graph.sub)$color
	# )
# }


## get OTUs into modules
membership <- membership(sparcc.comm.between.weighted)
modules <- ldply(membership, .id="OTU")
modules.sort <- arrange(modules, V1)
modules.sort

# rename module OTUs back to contain leading zeros (original format)
names <- as.character(modules.sort$OTU)
sample.split <- strsplit(names,'Otu',fixed=TRUE)
sample.split.temp <- ldply(sample.split)
sample.split.temp$V2 <- sprintf("%06d", as.numeric(sample.split.temp$V2))
names.num <- sample.split.temp$V2
names <- paste("Otu", names.num, sep="")

modules.sort$OTU <- names

otu.sub <- data.frame(OTU=rownames(tax_table(y.sum)), tax_table(y.sum))
str(otu.sub)
otu.sub.mod <- join(otu.sub, modules.sort, by="OTU")
otu.sub.mod.sort <- arrange(otu.sub.mod, V1)

dynamicColors <- labels2colors(modules$V1)
table(dynamicColors)
str(dynamicColors)


# print resulting graph with communities
sizeGrWindow(7, 7.77)

colbar <- names(table(dynamicColors))

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}



#################################
### OUTPUT ###
#################################
l.f <- layout.auto(sparcc.graph.sub)

# unconnected graph, vertices only
pdf("~/Dropbox/Work/grad_school/1_Baikal/5_Data Analysis/figures/network_graph.pdf", width=7, height=7.77)
plot.igraph(sparcc.graph.sub,
	layout=l.f,
	edge.color="white",
	vertex.size=exp(range01(V(sparcc.graph.sub)$weight))*3,
	vertex.color="black",
	vertex.label=NA)
dev.off()

# connected graph, vertices and edges
pdf("~/Dropbox/Work/grad_school/1_Baikal/5_Data Analysis/figures/network_graph_edge.pdf", width=7, height=7.77)
plot.igraph(sparcc.graph.sub,
	layout=l.f,
	edge.color=E(sparcc.graph.sub)$color,
	vertex.size=exp(range01(V(sparcc.graph.sub)$weight))*3,
	vertex.color="black",
	vertex.label=NA)
dev.off()

# connected graph with communities, vertices and edges
pdf("~/Dropbox/Work/1_baikal/5_data_analysis/figures/network_graph_edge_comm.pdf", width=7, height=7.77)
plot.communities(sparcc.comm.between.weighted, sparcc.graph.sub,
#		main="sparcc.comm.between.weighted",
		layout=l.f,
		vertex.size=exp(range01(V(sparcc.graph.sub)$weight))*3,
		# vertex.label=V(sparcc.graph.sub)$name,
		# vertex.label.cex=0.5,
#		vertex.label=NA,
		colbar=colbar,
		mark.col=add.alpha(col=colbar, alpha=0.3),
		mark.border=add.alpha(col=colbar, alpha=0.5),
		edge.color=E(sparcc.graph.sub)$color
)
dev.off()


### EIGENOTUS ###
MEList = moduleEigengenes(subset(t(otu_table(y.sum)), select=modules.sort$OTU), colors = dynamicColors)
MEs = MEList$eigengenes
colnames(MEs)

modules <- data.frame(MEs, sample_data(y.sum))
modules.plot <- melt(modules, measure.vars=colnames(MEs), variable.name="module", value.name="eigen")

MElist <- colnames(MEs)
#MElist <- c("MEblack", "MEturquoise", "MEred", "MEgreen", "MEblue", "MEbrown", "MEyellow")
mod.plot <- list()
MEgrad <- str_sub(MElist, 3, 100)

g1 <- ggplot(data=modules.plot, aes(x=temp, y=eigen, color=module)) +
	theme_bw() +
	scale_colour_manual(values=MEgrad, "Module\nEigenvector") +
	xlab(expression(paste("Temperature ",degree,"C"))) +
	geom_point() +
	stat_smooth(method="loess", se=F) +
	facet_wrap(~module)

ggsave(g1, file="~/Dropbox/Work/grad_school/1_Baikal/5_Data Analysis/figures/modules_temp.pdf")

modules.plot.vars <- melt(modules.plot, measure.vars=c("TN", "DN", "TP", "DP", "DS", "temp"), variable.name="vars", value.name="var.value")

ggplot(data=subset(modules.plot.vars, module=="MEgreen" & region!="PB"), aes(x=var.value, y=eigen, color=vars)) +
	geom_point() +
	stat_smooth(method="lm", se=F) +
	facet_wrap(~vars, scales="free_x")
	
	
	
	
	
	
	
	
	
	
	
	
	