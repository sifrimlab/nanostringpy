if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes most up to date version of Bioc
BiocManager::install(version="3.15")

BiocManager::install("NanoStringNCTools")
BiocManager::install("GeomxTools")
BiocManager::install("GeoMxWorkflows")


library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)


if(packageVersion("GeomxTools") < "2.1" & 
   packageVersion("GeoMxWorkflows") >= "1.0.1"){
  stop("GeomxTools and Workflow versions do not match. Please use the same version. 
    This workflow is meant to be used with most current version of packages. 
    If you are using an older version of Bioconductor please reinstall GeoMxWorkflows and use vignette(GeoMxWorkflows) instead")
}

if(packageVersion("GeomxTools") > "2.1" & 
   packageVersion("GeoMxWorkflows") <= "1.0.1"){
  stop("GeomxTools and Workflow versions do not match. 
         Please use the same version, see install instructions above.")
  
  # to remove current package version
  # remove.packages("GeomxTools")
  # remove.packages("GeoMxWorkflows")
  # see install instructions above 
}

datadir <- file.path("/home/david/Documents/als_overlay/nanostring/M-558KULeuven-Cruz/Raw_data_files/WTA/")
DCCFiles <- dir(file.path(datadir, "DCC"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)

PKCFiles <- dir(file.path(datadir, "Mm_R_NGS_WTA_v1.0/"), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)


SampleAnnotationFile <- dir(file.path(datadir), pattern = "template.xlsx$",full.names = TRUE, recursive = TRUE)

demoData <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "Template",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"),
                         experimentDataColNames = c("panel"))
demoData
library(knitr)
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

library(dplyr)
library(ggforce)

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols

count_mat <- count(pData(demoData), `slide name`, class,  region, segment)
nrow(pData(demoData))
count_mat
# gather the data and plot in order: class, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
test_gr
test_gr$x <- factor(test_gr$x,
                    levels = c("class", "region" , "segment", "slide name"))
# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = class), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") 

#########################
## QC & Pre-processing ##
#########################
# Shift all values to 0-1
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

## Segment QC
# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))


library(ggplot2)

#col_by <- "segment"
col_by <- "region"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}
QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
# QC_histogram(sData(demoData), "nuclei", col_by, 20) # not present in dataset

# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans
# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
col_by = "region"
for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}
# remove neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

sData(demoData)
# show all NTC values, Freq = # of Segments with a given NTC count:
# NTC is the no template control, whatever that is
kable(table(NTC_Count = sData(demoData)$NTC),
      col.names = c("NTC Count", "# of Segments"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")
demoData <- demoData[, QCResults$QCStatus == "PASS"]
dim(demoData)
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)
dim(demoData)
ProbeQCResults <- fData(demoData)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df

ProbeQCPassed <- 
  subset(demoData, 
         fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

demoData <- ProbeQCPassed 

length(unique(featureData(demoData)[["TargetName"]]))
target_demoData <- aggregateCounts(demoData)
dim(target_demoData)
exprs(target_demoData)[1:5, 1:2]

## QC based on LOQ
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ

# Filter out
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]
fData(target_demoData)$TargetName
sum(LOQ_Mat == 1)
sum(LOQ_Mat == 0)

LOQ_Mat

# pData(target_demoData)$LOQ
# Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = region)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

#########################
## Problem: selecting on 0.1 detection rate 
## Removes all "accidentally" removes all SN samples
##########################

kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$class))

dim(target_demoData)
target_demoData <-
  target_demoData[, pData(target_demoData)$GeneDetectionRate >= .1]

dim(target_demoData)


# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols

count_mat <- count(pData(target_demoData), `slide name`, class, region, segment)
# gather the data and plot in order: class, region, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
test_gr$x <- factor(test_gr$x,
                    levels = c("class",  "region" , "segment", "slide name"))
# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = class), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "")# +
  #annotate(geom = "segment", x = 4.25, xend = 4.25, y = 20, 
  #         yend = 120, lwd = 2) #+
  #annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
   #        hjust = 0.5, label = "100 segments")


# Gene Detection rate
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))
"Chat" %in% rownames(fData(target_demoData))
# Gene of interest detection table
goi <- c("Chat","Slc18a3","Ahnak2","Chodl","Slc5a7","Mmp9","Pdgfd","Grin3b","Isl1","Mnx1","Vipr2", "Fus", "Tardbp")

goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_demoData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_demoData)[goi, "DetectionRate"]))

goi_df

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_demoData <- 
  target_demoData[fData(target_demoData)$DetectionRate >= 0.1 |
                    fData(target_demoData)$TargetName %in% neg_probes, ]
dim(target_demoData)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_demoData)]
goi


###################
## Normalization ##
###################
library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "segment"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_demoData)),
             Segment = colnames(exprs(target_demoData)),
             Annotation = pData(target_demoData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_demoData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_demoData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_demoData <- normalize(target_demoData ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_demoData <- normalize(target_demoData ,
                             norm_method = "neg", 
                             fromElt = "exprs",
                             toElt = "neg_norm")

# visualize the first 10 segments with each normalization method
boxplot(exprs(target_demoData)[,1:12],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:12, xlab = "Segment",
        ylab = "Counts, Raw")
boxplot(assayDataElement(target_demoData[,1:12], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:12, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

boxplot(assayDataElement(target_demoData[,1:12], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:12, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")


#################
## Unsupervised #
#################
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),  
       config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = class, shape = segment)) +
  geom_point(size = 3) +
  theme_bw()


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_demoData , elt = "q_norm"))),
        perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = class, shape = segment)) +
  geom_point(size = 3) +
  theme_bw()

########### Clustering high CV genes
library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_demoData, elt = "log_q") <-
  assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_demoData,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
pheatmap(assayDataElement(target_demoData[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           #pData(target_demoData)[, c("class", "segment","region")])
           pData(target_demoData)[, c("class", "segment")])
           #pData(target_demoData)[, c( "segment")])

#################
## DE analysis ##
#################
# convert test variables to factors
pData(target_demoData)$testRegion <- 
  #factor(pData(target_demoData)$class, c("CTRL", "FUS", "TDP-43"))
  factor(pData(target_demoData)$segment, c("ChATpos", "ChATneg"))
pData(target_demoData)[["slide"]] <- 
  factor(pData(target_demoData)[["replicate_nr"]])
assayDataElement(object = target_demoData, elt = "log_q") <-
  assayDataApply(target_demoData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
#pData(target_demoData[,ind])[["slide"]]
#length(pData(target_demoData[,ind]))
#pData(target_demoData[,ind])
for(status in c("CTRL", "FUS", "TDP-43")) {
  ind <- pData(target_demoData)$class == status
  mixedOutmc <-
    mixedModelDE(target_demoData[, ind],
                 #target_demoData,
                 elt = "log_q",
                 #modelFormula = ~ testRegion + (1 + testRegion | slide),
                 modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion",
                 #nCores = parallel::detectCores())
                 nCores = 1,
                 multiCore = TRUE)
  
   #format results as data.frame
   r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
   tests <- rownames(r_test)
   r_test <- as.data.frame(r_test)
   r_test$Contrast <- tests
   
   # use lapply in case you have multiple levels of your test factor to
   # correctly associate gene name with it's row in the results table
   r_test$Gene <- 
     unlist(lapply(colnames(mixedOutmc),
                   rep, nrow(mixedOutmc["lsmeans", ][[1]])))
   r_test$Subset <- status
   r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
   r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                        "Pr(>|t|)", "FDR")]
   results <- rbind(results, r_test)
}

### interpreting results table
kable(subset(results, Gene %in% goi), digits = 3,
      caption = "DE results for Genes of Interest",
      align = "lc", row.names = FALSE)

###############################################
## Between Slide Analysis: Disease vs Healty ##
###############################################
# convert test variables to factors
pData(target_demoData)$testClass <-
  factor(pData(target_demoData)$class, c("CTRL", "FUS", "TDP-43"))

# run LMM:
# formula follows conventions defined by the lme4 package
results2 <- c()
for(segment in c("ChATpos", "ChATneg")) {
  ind <- pData(target_demoData)$segment == segment
  mixedOutmc <-
    mixedModelDE(target_demoData,
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 | slide),
                 groupVar = "testClass",
                 nCores = 1,
                 multiCore = TRUE)
                 #nCores = parallel::detectCores(),
                 #multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test
  r_test$Subset <- segment
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                     "Pr(>|t|)", "FDR")]
  r_test
  results2 <- rbind(results2, r_test)
}
results2
    
kable(subset(results2, Gene %in% goi), digits = 3,
    caption = "DE results for Genes of Interest",
    align = "lc", row.names = FALSE)

###### Visualizing DE gene
library(ggrepel) 
# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                      levels = c("NS or FC < 0.5", "P < 0.05",
                                 "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
#for(cond in c("ChATneg", "ChATpos")) {
for(cond in c("CTRL", "FUS", "TDP-43")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
               results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix
# 
length(results$Gene)
length(unique(results$Gene))
#tmp_top_g <- top_g[0:3]
# Graph results
ggplot(results,
     aes(x = Estimate, y = -log10(`Pr(>|t|)`),
         color = Color, label = Gene)) +
geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
geom_hline(yintercept = -log10(0.05), lty = "dashed") +
geom_point() +
labs(x = "log2(FC)",
     y = "Significance, -log10(P)",
     color = "Significance") +
scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                              `FDR < 0.05` = "lightblue",
                              `P < 0.05` = "orange2",
                              `NS or FC < 0.5` = "gray"),
                   guide = guide_legend(override.aes = list(size = 4))) +
scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
geom_text_repel(data = subset(results, Gene %in% top_g), #& FDR < 0.001),
                size = 4, point.padding = 0.15, color = "black",
                min.segment.length = .1, box.padding = .2, lwd = 2,
                max.overlaps = 50) +
theme_bw(base_size = 16) +
theme(legend.position = "bottom") +
facet_wrap(~Subset, scales = "free_y")


### Plot Genes of interest
kable(subset(results, Gene %in% goi), row.names = FALSE)

ggplot(pData(target_demoData),
     aes(x = class, fill = segment,
         y = assayDataElement(target_demoData["Chat", ],
                              elt = "q_norm"))) +
geom_violin() +
geom_jitter(width = .2) +
labs(y = "Chat Expression") +
scale_y_continuous(trans = "log2") +
facet_wrap(~region) +
theme_bw()

####################################
## Plot targets against eachother ##
####################################

glom <- pData(target_demoData)$region == "ChATpos"

gene1 <- "Fus"
gene2 <- "Tardbp"
# show expression of PDHA1 vs ITGB1
ggplot(pData(target_demoData),
       aes(x = assayDataElement(target_demoData[gene1, ],
                                elt = "q_norm"),
           y = assayDataElement(target_demoData[gene2, ],
                                elt = "q_norm"),
           color = region)) +
  geom_vline(xintercept =
               max(assayDataElement(target_demoData[gene1, glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_hline(yintercept =
               max(assayDataElement(target_demoData[gene2, !glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  labs(x = sprintf("%s Expression", gene1), y = sprintf("%s Expression", gene2)) +
  facet_wrap(~class)



##################################
## Heatmap of significant genes ##
##################################

#select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(results, `Pr(>|t|)` < 0.05)$Gene)
GOI
pheatmap(log2(assayDataElement(target_demoData[GOI, ], elt = "q_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 2, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_demoData)[, c("segment", "class")])
