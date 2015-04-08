##
#Network analysis with igraph
#tab-delimited text file for input network_fp
#must have columns: "var1", "var2", "weight", "pvalue"

igraph.f=function(network_fp,name, pv_threshold=0.050){
  library(igraph)
  data=read.table(network_fp, header=TRUE, check.names=FALSE, sep="\t")
  data2=data[(2*is.na(data[,"pvalue"])==0),]
  data2.5=data2[round(data2[,"pvalue"],digits=3)<=pv_threshold,]
  data3=data2.5[,1:3]
  colnames(data3)=c("var1", "var2", "weight")
  g=graph.data.frame(data3,directed=FALSE)
  
  
  #diameter
  d=diameter(g)
  
  print(summary.igraph(g))
  #no. nodes
  n=g[1]
  
  #no edges
  m=length(g[3][[1]])
  
  
  #number of cliques?
  cliq=length(cliques(g))
  
  #mean closensss (average length of shortest path to/from all other vertices in graph)
  #small world effect, l
  l=(closeness(g, v=V(g), mode="all"))
  l.mean=mean(l)
  l.max=max(l)
  l.min=min(l)
  #
  
  #transitivity, C (same as clustering coefficient in cytoscape)
  #C is average
  C=transitivity(g, type="global")
  #Ci is a list for every vertex
  Ci=transitivity(g, type="local")
  
  #degree distribution
  pk=degree(g, v=V(g), mode="all")
  pk.max=max(pk)
  pk.min=min(pk)
  z=mean(pk)
  
  
  #power law exponent, alpha, fit to the degree distribution
  p=power.law.fit(pk)
  alpha=p$coef
  
  #power law exponent, alpha, fit to the cumulative degree distribution, so add 1
  dd=degree.distribution(g, cumulative=TRUE)
  power.law.fit(dd+1)
  
  #degree correlation coefficient: (do high/low -degree nodes preferentially associated with each other?, Pearson's correlation approach where positive number is assortative, negative is disassortative, from Table 2, value   in Newman
  node.names=as.vector(unlist(g[[9]][[3]]))
  node.degrees=pk
  edgelist=get.edgelist(g)
  
  degree.edge1=NULL
  for(i in 1:length(edgelist[,1])){
    degree.edge1[i]=node.degrees[edgelist[i,1]==node.names]
  }
  
  degree.edge2=NULL
  for(i in 1:length(edgelist[,2])){
    degree.edge2[i]=node.degrees[edgelist[i,2]==node.names]
  }
  
  PearsCor=cor.test(degree.edge1, degree.edge2, method="pearson", conf.level=0.95)
  r=PearsCor$estimate
  r.p=PearsCor$p.value
  
  
  #another possible function that calculates similarity
  test=similarity.dice(g, vids=V(g), mode="all")
  test.d=as.dist(test)
  
  
  #insert assortivity (mixing) by phylum-level assignment, equ.17 in Newman
  
  pdf(paste("MINE_",name,"_igraphPlots.pdf",sep=""),onefile=TRUE)
  par(mar=c(6,4,2,1)+0.1)
  par(mfrow=c(2,2))
  plot(degree.distribution(g, cumulative=TRUE), main="Degree Distributions", xlab="No. edges", ylab="Cumulative proportion of nodes", type="b", col="red")
  hist(l, breaks=20, col="gray", xlab="Mean geodesic distance")
  plot(degree.edge2,degree.edge1, main="Degree correlation")
  plot(test.d, main="Dice similarity between all pairs of Node")
  dev.off()
  
  #output
  out=c(n,m,z,d,l.mean,alpha,C,r, r.p)
  out=unlist(out)
  names(out)=c("No_Nodes_n", "No_Edges_m","Mean_degree","Diameter", "Mean_geodesic_l","PwrLawAlpha", "ClusteringCoef_C", "PearsonR", "Pearson_pvalue")
  
  write.table(t(out), paste("MINE_",name, "_MIC_networkStats.txt",sep=""), sep="\t", quote=FALSE, row.names=FALSE)
  print(out)
  return(out)
}

#Load the above function, then use
#network_fp="TYPE_FILE_NAME"
igraph.f(network_fp,name, pv_threshold)
