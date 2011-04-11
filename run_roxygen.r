library(roxygen)

setwd("~/Raivo/Projects/RHeatmap/")

roxygenize(package.dir = "Source", roxygen.dir = "pheatmap", unlink.target = T, use.Rd2 = T)

## Run in terminal
cd ~/Raivo/Projects/RHeatmap/
R CMD check pheatmap
R CMD build pheatmap
R CMD install pheatmap_0.4.tar.gz 


# Debug, kui kisab et dokument sisaldab mitte ascii tähti
tools::showNonASCII( readLines("RUtil/man/gprofiler.Rd"))