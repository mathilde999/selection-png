#!/usr/bin/Rscript

#Intrapolate recombination rate from reference (1000g or hapmap files)

usage <- paste('usage: Rscript GeneticMapIntrapolate.r <plinkFile> <outFile> <ref_map>\n')
if (!exists("parameters"))
{
    parameters <- commandArgs(trailingOnly=T)
}
if (length(parameters) == 3)
{
    plinkFile <- parameters[1]
    outFile <- parameters[2]
    ref_map <- parameters[3]
} else
{
    stop(usage)
}
mapFile <- read.table(file = plinkFile, header = FALSE)
names(mapFile) <- c('chr', 'snpId', 'gd', 'pos')
gd <- read.table(file = ref_map, sep = " ")
names(gd) <- c('chr', 'snpId', 'gd', 'pos')
hapMapPositions <- gd$pos
minGd <- min(hapMapPositions)
maxGd <- max(hapMapPositions)
genMaps<- gd$gd
previous_positions <- mapFile$pos[which(mapFile$pos <= minGd)]
posterior_positions <- mapFile$pos[which(mapFile$pos >= maxGd)]
interpolating_positions <-  mapFile$pos[intersect(which(mapFile$pos > minGd), which(mapFile$pos < maxGd))]
genMapsInterpolated <-approx(x = hapMapPositions, y = genMaps, xout = interpolating_positions, rule = 2)$y
mapFile$gd[ which(mapFile$pos <= minGd)] <- 0
print(head(mapFile))
mapFile$gd[which(mapFile$pos >= maxGd)] <- max(genMaps)
mapFile$gd[ intersect(which(mapFile$pos > minGd), which(mapFile$pos < maxGd))] <- genMapsInterpolated
write.table(mapFile, outFile, quote = FALSE, row.names = FALSE,col.names=FALSE,sep="\t")

