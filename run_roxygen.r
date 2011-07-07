library(roxygen)

setwd("~/Raivo/Projects/RHeatmap/")

roxygenize(package.dir = "Source", roxygen.dir = "pheatmap", unlink.target = T, use.Rd2 = T)

## Run in terminal
cd ~/Raivo/Projects/RHeatmap/
rm -r pheatmap/.git
rm -r pheatmap/inst
R CMD check pheatmap
R CMD build pheatmap
R CMD install pheatmap_0.5.2.tar.gz 


# Debug, kui kisab et dokument sisaldab mitte ascii tähti
tools::showNonASCII( readLines("RUtil/man/gprofiler.Rd"))