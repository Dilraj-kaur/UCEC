source("http://bioconductor.org/biocLite.R")
biocLite("hgu133a.db")
library(hgu133a.db)
library(annotate)
x <- hgu133aSYMBOL
# Get the probe identifiers - gene symbol mappings
mapped_probes <- mappedkeys(x)
# Convert to a dataframe
genesym.probeid <- as.data.frame(x[mapped_probes])
head(genesym.probeid)