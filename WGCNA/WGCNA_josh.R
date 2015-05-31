### WGCNA

[WGCNA Home Page](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/)

[WGCNA R Tutorial](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf)

### R code for WGCNA

# this is the R code I used to install WGCNA

getCRANmirrors(all = FALSE, local.only = FALSE)
options(repos=structure(c(CRAN="http://ftp.ussg.iu.edu/CRAN/")))

source("http://bioconductor.org/biocLite.R")

biocLite("bioDist")
biocLite("impute")
biocLite("preprocessCore")
install.packages("WGCNA")

# analysis code is a WIP