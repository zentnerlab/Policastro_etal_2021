library(ZebrafishDevelopmentalCAGE)
library(CAGEr)
library(writexl)
library(GenomicFeatures)
library(GenomicRanges)
library(ChIPseeker)

# Load data
data(ZebrafishSamples)
data(ZebrafishCAGE)

# Import samples 
zf_cageset <- importPublicData(source = "ZebrafishDevelopment", 
                               dataset = "ZebrafishCAGE", 
                               group = "development", 
                               sample = c("zf_unfertilized_egg", "zf_prim20"))

# Normalize
plotReverseCumulatives(zf_cageset, fitInRange = c(5,1000), onePlot = TRUE)
dev.off()

normalizeTagCount(zf_cageset, method = "powerLaw", fitInRange = c(5,1000), alpha = 1.18, T = 1e6)

# Cluster TSSs into TSRs
clusterCTSS(object = zf_cageset, threshold = 3, thresholdIsTpm = TRUE, nrPassThreshold = 1,
            method = "distclu", maxDist = 25, removeSingletons = TRUE, keepSingletonsAbove = 3,
            useMulticore = TRUE)

# Cluster TSRs into a consensus set
aggregateTagClusters(zf_cageset, tpmThreshold = 10, maxDist = 100)

# Get cumulative CAGE distributions over consensus TSRs
cumulativeCTSSdistribution(zf_cageset, clusters = "consensusClusters", useMulticore = TRUE)

# Calculate shifting scores
scoreShift(zf_cageset, groupX = "zf_unfertilized_egg", groupY = "zf_prim20", testKS = TRUE, useTpmKS = FALSE)

# Retrieve shifts without score thresholding
shifting.promoters.noScoreThresh <- getShiftingPromoters(zf_cageset,
                                                         tpmThreshold = 10,
                                                         scoreThreshold = -Inf,
                                                         fdrThreshold = 0.05)

# Retrieve shifts with score thresholding
shifting.promoters.scoreThresh.0.4 <- getShiftingPromoters(zf_cageset,
                                                           tpmThreshold = 10,
                                                           scoreThreshold = 0.4,
                                                           fdrThreshold = 0.05)

# Annotate shifted TSRs
annotation <- makeTxDbFromGFF("danRer7.refGene.gtf")

no_thresh_gr <- makeGRangesFromDataFrame(shifting.promoters.noScoreThresh, keep.extra.columns = TRUE)
thresh_04_gr <- makeGRangesFromDataFrame(shifting.promoters.scoreThresh.0.4, keep.extra.columns = TRUE)

no_thresh_anno <- annotatePeak(no_thresh_gr, TxDb = annotation, tssRegion = c(-500,500), sameStrand = TRUE)
thresh_04_anno <- annotatePeak(thresh_04_gr, TxDb = annotation, tssRegion = c(-500,500), sameStrand = TRUE)

# Export shifts to XLSX
write_xlsx(list("No score threshold" = as.data.frame(no_thresh_anno@anno), 
                "Score threshold = 0.4" = as.data.frame(thresh_04_anno@anno)),
           path = "CAGEr_shifts.xlsx")
