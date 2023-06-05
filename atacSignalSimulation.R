#@' author Emanuel Sonder
library(Rsamtools)
library(GenomicAlignments)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MASS)
library(mclust)
library(Repitools)

first <- data.table::first

sampleSwitch <- function(total, size){
  if(total >= size)
  {
    s <- sample(total, size=size, replace=FALSE)
  }
  else
  {
    s <- sample(total, size=size, replace=TRUE)
  }
  
  return(s)
}

.importFrags <- function(bamPath, which=NULL, annotationStyle="NCBI")
{
  if(is.null(which))
  {
    which <- GRanges(Rle(1:22), 
                     IRanges(start=rep(1,22), width=536870912))
  }
  
  param <- ScanBamParam(what=c('pos', 'qwidth', 'isize'),
                        which=which, 
                        flag=scanBamFlag(isUnmappedQuery=FALSE))
  
  readPairs <- readGAlignmentPairs(bamPath, param=param)
  
  # get fragment coordinates from read pairs
  frags <- GRanges(seqnames(GenomicAlignments::first(readPairs)), 
                   IRanges(start=pmin(GenomicAlignments::start(GenomicAlignments::first(readPairs)), 
                                      GenomicAlignments::start(GenomicAlignments::second(readPairs))), 
                           end=pmax(GenomicAlignments::end(GenomicAlignments::first(readPairs)), 
                                    GenomicAlignments::end(second(readPairs)))))
  
  frags <- granges(frags, use.mcols=TRUE)
  
  # change annotation style
  seqlevelsStyle(frags) <- annotationStyle
  #print(length(frags))
  
  return(list(fragments=frags,
              readPairs=readPairs))
}

.importPeaks <- function(peaksDir, 
                         which=NULL, 
                         annotationStyle="NCBI")
{
  peaks <- fread(peaksDir, select=1:3)
  colnames(peaks) <- c("chr", "start", "end")
  
  peaks <- makeGRangesFromDataFrame(peaks)
  seqlevelsStyle(peaks) <- annotationStyle
  
  peaks <- subsetByOverlaps(peaks, which)
  
  return(peaks)
}


# estimate the log fold change based on the enrichment over the input
.estLfc <- function(enrichment, 
                    lfcDist,
                    nBins=100, 
                    enrCol="rel_enr")
{
  enrPeakDt <- as.data.table(enrichment)
  lfcDistDt <- as.data.table(lfcDist)
  
  # bin enrichment values
  bins <- seq(min(lfcDistDt[[enrCol]]), 
              max(lfcDistDt[[enrCol]]), length.out=nBins)
  lfcDistDt[,bin:=findInterval(get(enrCol), bins)]
  lfcDistDt[,bin_start:=bins[bin]]
  lfcDistDt[,bin_end:=bins[bin+1]]
  lfcDistDt[is.na(bin_end),]$bin_end <- 100
  lfcDistDt[bin_start==min(bin_start),]$bin_start <- -100
  
  # sample values from respective bins
  # overlap by which bin contains 
  enrPeakDt[,peak_id:=1:nrow(enrPeakDt)]
  enrPeakDt[,start:=get(enrCol)]
  enrPeakDt[,end:=get(enrCol)]
  
  setkey(lfcDistDt, bin_start, bin_end)
  enrPeakDt <- foverlaps(enrPeakDt, 
                         lfcDistDt[,c("bin_start", "bin_end", "lfc")], 
                         by.x=c("start", "end"),
                         by.y=c("bin_start", "bin_end"),
                         type=c("any"))
  enrPeakDt$start <- NULL
  enrPeakDt$end <- NULL
  
  # sample one value per bin
  enrPeakDt <- enrPeakDt[,.SD[sample(.N, 1)],by=peak_id]
  enrPeakDt[,lfc:=fifelse(is.na(lfc), 0, lfc)]
  setorder(enrPeakDt, peak_id)
  
  return(enrPeakDt)
}

.varSignal <- function(frags, 
                       peaks, 
                       frac=0.2, 
                       minOverlap=1, 
                       annotationStyle="NCBI")
{
  # Subset to "primary" signal peaks
  fragsPeaks <- suppressWarnings(subsetByOverlaps(frags, peaks, minoverlap=minOverlap))
  fragsOut <- suppressWarnings(subsetByOverlaps(frags, peaks, invert = TRUE))
  
  # Subsample reads outside of peaks
  nOut <- round((length(fragsPeaks)/(1-frac))*frac)
  fragsOutSub <- fragsOut[sample(1:length(fragsOut), min(nOut, length(fragsOut)))]
  
  fragsSub <- c(fragsPeaks, fragsOutSub)
  
  return(fragsSub)
}

.varyGCBias <- function(frags, 
                        biasFileDir, 
                        frac=0.8, 
                        minGC=0, 
                        maxGC=1,
                        annotationStyle="NCBI",
                        genome=BSgenome.Hsapiens.UCSC.hg38)
{
  # Import GC bias table (output deeptools)
  gcBiasTable <- fread(biasFileDir, select=1:3,
                       col.names=c("observed_read_count", 
                                   "expected_read_count", 
                                   "ratio"))
  gcBiasTable[,gc_bin:=seq(from=0,to=1, length.out=nrow(gcBiasTable))]
  
  # subset gc rates
  gcBiasTableSub <- subset(gcBiasTable, gc_bin>=minGC & gc_bin<=maxGC)
  
  gcBiasTableSub[,gc_bin_start:=gc_bin]
  gcBiasTableSub[,gc_bin_end:=data.table::shift(gc_bin_start, 1, type="lead")]
  gcBiasTableSub[nrow(gcBiasTableSub),]$gc_bin_end <- maxGC
  
  # get gc content per frag
  seqlevelsStyle(frags) <- "UCSC"
  frags$GC_content <- suppressWarnings(Repitools::gcContentCalc(frags, 
                                                                organism=genome))
  seqlevelsStyle(frags) <- annotationStyle
  
  # annotate each frag with gc bin
  fragsTable <- as.data.table(frags) # convert to data.table for simplicity
  fragsTable$min_GC_content <- fragsTable$GC_content
  fragsTable$max_GC_content <- fragsTable$GC_content  
  
  setkey(gcBiasTableSub, gc_bin_start, gc_bin_end)
  fragsTable <- foverlaps(fragsTable, gcBiasTableSub,
                          by.x=c("min_GC_content", "max_GC_content"),
                          by.y=c("gc_bin_start", "gc_bin_end"))
  
  # Calculate probability for sampling for each gc bin
  nSample <- round(nrow(fragsTable)*frac)
  gcBinObserved <- fragsTable[,.(ratio=unique(ratio),
                                 expected_read_count=unique(expected_read_count)), 
                              by=gc_bin]
  gcBinObserved[,norm_observed:=expected_read_count*ratio]
  gcBinObserved[,prob:=norm_observed/sum(norm_observed)]
  gcBinObserved[, n2Sample:=round(prob*nSample)]
  
  fragsTable <- merge(fragsTable, gcBinObserved, by=c("gc_bin"))
  
  # sample from gc bins
  # alternative with no up-sampling : sample(.N, min(unique(n2Sample), .N))
  fragsSubTable <- fragsTable[,.SD[sampleSwitch(.N, unique(n2Sample))], by=gc_bin]
  
  return(fragsSubTable)
}

.varyFragDist <- function(fragsTable, 
                          frac=0.8,
                          nClust=4, 
                          fitGMM=FALSE,
                          referenceData=NULL,
                          prob=c(0.5, 0.4, 0.08, 0.02),
                          estimateProb=FALSE)
{
  if(fitGMM)
  {
    # Estimate fragment length distributions
    if(is.null(referenceData))
    {
      fit = Mclust(fragsTable$width, G=nClust, model="V")
    }
    else
    {
      fit = Mclust(referenceData$width, G=nClust, model="V")
    }
    
    # probensities should be able to be changed by the users
    if(estimateProb)
    {
      prob <- fit$parameters$pro
    }
    
    # Annotate each fragment with the most likely cluster
    fragsTable$cluster <- fit$classification
  }
  else
  {
    fragsTable[, cluster:=fifelse(width<120, 1, as.integer(NA))]
    fragsTable[, cluster:=fifelse(width>=120 & width<300, 2, cluster)]
    fragsTable[, cluster:=fifelse(width>=300 & width<500, 3, cluster)]
    fragsTable[, cluster:=fifelse(width>=500, 4, cluster)]
  }
  
  # Number of reads to sample for each type
  nSample <- round(nrow(fragsTable)*frac)
  nFragsType <- round(nSample*prob)
  
  fragsSubTable <- lapply(1:nClust, function(i){
    sampledFrags <- fragsTable[cluster==i, ][sampleSwitch(.N, nFragsType[i]),]
  })
  
  # Exact but slow -------------------------------------------------------------
  # Annotate each fragment with the probability of originating from 
  # either NF, 1-, 2-, 3- nucleosome containing fragments
  # readsTable <- cbind(readsTable, fit$z)
  
  # Check if it has already been sampled
  # readsTable[,id:=1:nrow(readsTable)]
  # readsTable[,isSampled:=FALSE]
  
  
  # Sample each read with the respective probability
  #readsSubTable <- lapply(1:nClust, function(i){
  #  sampledReads <- readsTable[isSampled!=TRUE, ][sample(.N, nReadsType[i], 
  #                                                       prob=get(paste0("V", i))),]
  #  readsTable[,isSampled:=fifelse(isSampled | (id %in% sampledReads$id), 
  #                                 TRUE, FALSE)]
  #  return(sampledReads)
  #})
  
  #-----------------------------------------------------------------------------
  
  fragsSubTable <- rbindlist(fragsSubTable)
  #fragsSubTable <- unique(fragsSubTable, by=c("seqnames", "start", "end"))
  
  return(fragsSubTable)
}

.varEffectSize <- function(fragsTable, 
                           peaks,
                           logFCs=NULL,
                           effectStrength=0.5,
                           noiseLevel=0.1
                           #mode=c("scale", "plain")
){
  
  #mode <- mode[1]
  
  if(is.null(logFCs))
  {
    logFCs <- rep(0, length(peaks))
  }
  
  
  # convert to data.table for simplicity
  peaksTable <- as.data.table(peaks)
  peaksTable$id <- 1:nrow(peaksTable)
  peaksTable$logFC <- logFCs
  
  fragsTable[,frag_id:=paste(seqnames, start, end, sep="_")]
  
  # get frags within peaks
  setkey(peaksTable, seqnames, start, end)
  fragsInPeaks <- foverlaps(fragsTable, 
                            peaksTable,
                            type="any",
                            nomatch=NULL,
                            by.x=c("seqnames", "start", "end"),
                            by.y=c("seqnames", "start", "end"),
                            minoverlap=1,
                            mult="all")
  
  # Get frags outside of peaks
  fragsOutPeaks <- fragsTable[!unique(fragsInPeaks, by="frag_id"), on=.(frag_id)]
  
  # Calculate the number of fragments to sample per peak
  fragsInPeaks[,fc:=2^abs(logFC)]
  fragsInPeaks[, nFrags:=.N*round(data.table::first(fc)*effectStrength), by=c("id")]
  
  # Only sample from positive logFcs run for multiple samples
  if(nrow(fragsInPeaks[logFC<0,])>0)
  {
    fragsInPeaksSub <- fragsInPeaks[logFC<0,][,.SD[sampleSwitch(.N, data.table::first(nFrags))], by=id]
  }
  else
  {
    fragsInPeaksSub <- fragsInPeaks[logFC<0,]
  }
  
  fragsInPeaksSub <- rbind(fragsInPeaksSub, fragsInPeaks[logFC>=0,])
  fragsInPeaksSub <- fragsInPeaksSub[,c("seqnames", "i.start", "i.end"), with=FALSE]
  colnames(fragsInPeaksSub) <- c("seqnames", "start", "end")
  
  fragsSubset <- rbind(fragsInPeaksSub, 
                       fragsOutPeaks[,  c("seqnames", "start", "end"), with=FALSE])
  
  return(fragsSubset)
}

#' Function to vary different ATAC-seq properties
#'
#' Function to vary different Atac-seq related properties by subsampling fragments.
#' It is possible to vary GC-Bias, Fragment size distribution and effect strenght (of an experimental condition)
#' 
#' @param bamPath path to the .Bam file to subsample fragments from
#' @param bedPath path to the .Bed file with the peaks
#' @param sampleName name of the .Bam file to be analyized (subsampled .Bam gets saved under this name)
#' @param design a vector of 1 & -1 indicating to which of the two groups a
#' sample belongs to
#' @param logFCs vector with log fold changes for each peak. Needs to be of the same length 
#' as the peaks in the bed file. 
#' @param which GRanges object stating which regions should be selected
#' @param fracSub proportion of reads to subsample in each step 
#' (three steps in total, GC-bias based, fragment size based and effect strength based subsampling)
#' @param fragFragsOutside proportion of fragments outside of peaks
#' @param effectStrength Factor to which the fragments of differentially accessible peaks should be
#' sub- / upsampled. 
#' @param nFragTypes Number of fragment types (nucleosome-free, mono-, di-, tri-nucleosome, etc. containing fragments)
#' @param fitGMM TRUE/FALSE if GMM should be fitted for fragment length distribution or fixed thresholds should be used.
#' @param prob vector of probabilities of fragment types, needs same length as nFragTypes
#' @param estimateProb TRUE/FALSE if probabilities of fragment types should be estimated from data
#' @param minOverlap minimal overlap to consider a fragment being within a peak
#' @param minGC minimal gc content of fragments to be considered (to be removed)
#' @param maxGC maximal gc content of fragments to be considered (to be removed)
#' @param annotationStyle Either NCBI or UCSC
#' @param genome BSgenome object to be used. 
#' @return data.table with subsampled fragment coordinates & .bam file of these fragments saved on disk.
varyAtacSignal <- function(bamPath, 
                           bedPath, 
                           biasFileDir,
                           sampleName,
                           logFCs,
                           which=NULL, 
                           fracSub=0.8,
                           effectStrength=0.5,
                           nFragTypes=4,
                           fitGMM=FALSE,
                           prob=c(0.6, 0.3, 0.07, 0.03),
                           estimateProb=FALSE,
                           minOverlap=1,
                           minGC=0,
                           maxGC=1,
                           annotationStyle="NCBI",
                           genome=BSgenome.Hsapiens.UCSC.hg38)
{
  # Import fragments and peaks
  bamData <- .importFrags(bamPath, which, annotationStyle)
  frags <- bamData$fragments
  readPairs <- bamData$readPairs
  peaks <- .importPeaks(bedPath, which, annotationStyle)
  
  # Vary GC Bias
  fragsSubset <- .varyGCBias(frags, biasFileDir, fracSub,
                             minGC, maxGC, annotationStyle, genome)
  
  # Vary Frag Dist
  fragsSubset <- .varyFragDist(fragsSubset, fracSub, nClust=nFragTypes,
                               fitGMM=fitGMM, prob=prob, 
                               estimateProb=estimateProb)
  
  # Vary Effect size
  fragsSubset <- .varEffectSize(fragsSubset, peaks, effectStrength, logFCs=logFCs)
  
  # Get indices of read pairs to keep
  readPairsFrag <- data.table(seqnames=runValue(seqnames(GenomicAlignments::first(readPairs))), 
                              start=pmin(start(GenomicAlignments::first(readPairs)), 
                                         start(GenomicAlignments::second(readPairs))), 
                              end=pmax(end(GenomicAlignments::first(readPairs)), 
                                       end(GenomicAlignments::second(readPairs))),
                              id=1:length(readPairs))
  
  fragsSubset[,frag_id:=1:nrow(fragsSubset)]
  readPairsFragSub <- merge(fragsSubset,
                            readPairsFrag,
                            by.x=c("seqnames", "start", "end"),
                            by.y=c("seqnames", "start", "end"),
                            allow.cartesian=TRUE, # sample with replacement can cause this behavior
                            all.x=FALSE)
  readPairsFragSub <- unique(readPairsFragSub, by=c("frag_id"))
  readPairsSub <- readPairs[readPairsFragSub$id]
  
  # Convert to GRanges and get GC content
  fragsSubset <- makeGRangesFromDataFrame(fragsSubset)
  seqlevelsStyle(fragsSubset) <- "UCSC"
  fragsSubset$gc_content <- Repitools::gcContentCalc(fragsSubset, 
                                                     organism=genome)
  
  #seqlevelsStyle(fragsSubset) <- annotationStyle
  
  return(fragsSubset)
}

#' Function to simulate a two condition ATAC seq experiment
#'
#' Single samples are varied by group specific parameters in a two group setting.
#' Fragments are sub-samples based on a GC-Bias, Fragment size distribution and
#' effect strength which can be varied (see varyAtacSignal). 
#' 
#' @param bamPaths a vector of the paths of Bam files to subsample from. 
#' @param bedPath path to the .Bed file with the peaks
#' @param sampleNames a vector with the names of Bam files analyzed
#' @param design a vector of 1 & -1 indicating to which of the two groups a
#' sample belongs to
#' @param paramsGroup1 data.table with parameters for group 1 (see function: varyAtacSignal)
#' @param paramsGroup2 data.table with parameters for group 2 (see function: varyAtacSignal)
#' @param logFCs vector with log fold changes for each peak. Needs to be of the same length 
#' as the peaks in the bed file. 
#' @param which GRanges object stating which regions should be selected
#' @param effectStrength Factor to which the fragments of differentially accessible peaks should be
#' sub- / upsampled. 
#' @param nFragTypes Number of fragment types (nucleosome-free, mono-, di-, tri-nucleosome, etc. containing fragments)
#' @param minOverlap minimal overlap to consider a fragment being within a peak
#' @param minGC minimal gc content of fragments to be considered (to be removed)
#' @param maxGC maximal gc content of fragments to be considered (to be removed)
#' @param annotationStyle Either NCBI or UCSC
#' @param genome BSgenome object to be used. 
#' @return data.table with subsampled fragment coordinates & .bam file of these fragments saved on disk.
simAtacData <- function(bamPaths, 
                        bedPath, 
                        sampleNames,
                        gcBiases, 
                        design,
                        paramsGroup1,
                        paramsGroup2, 
                        enrichment,
                        lfcDist,
                        which=NULL,
                        effectStrength=2,
                        nFragTypes=4,
                        minOverlap=1,
                        minGC=0,
                        maxGC=1,
                        annotationStyle="NCBI",
                        genome=BSgenome.Hsapiens.UCSC.hg38){
  
  # seq level style
  seqlevelsStyle(enrichment) <- annotationStyle
  enrichment <- subsetByOverlaps(enrichment, which)
  
  # estimate logFCs
  peakDt <- .estLfc(enrichment, lfcDist)
  logFCs <- peakDt$lfc
  print(length(unique(logFCs)))
  
  # Positive samples
  posSamples <- bamPaths[which(design==1)]
  posSampleNames <- sampleNames[which(design==1)]
  posGcBiases <- gcBiases[which(design==1)]
  
  simSamplesPos <- lapply(1:length(posSamples), function(i){
    
    simData <- varyAtacSignal(bamPath=posSamples[i], 
                              bedPath=bedPath, 
                              biasFileDir=posGcBiases[i],
                              sampleName=posSampleNames[i],
                              logFCs=logFCs,
                              effectStrength=effectStrength,
                              nFragTypes=nFragTypes,
                              fracSub=paramsGroup1$fracSub[1],
                              prob=c(paramsGroup1$prob_nf[1], 
                                     paramsGroup1$prob_mono[1], 
                                     paramsGroup1$prob_di[1], 
                                     paramsGroup1$prob_tri[1]),
                              estimateProb=paramsGroup1$estimateProb,
                              which=which,
                              minOverlap=minOverlap,
                              minGC=minGC,
                              maxGC=maxGC,
                              annotationStyle=annotationStyle,
                              genome=genome)
    simData$sample <- posSampleNames[i]
    simData$group <- paramsGroup1$name[1]
    simData <- as.data.table(simData)
    
    return(simData)
  })
  
  simSamplesPos <- rbindlist(simSamplesPos)
  
  # Negative samples
  negSamples <- bamPaths[which(design==-1)]
  negSampleNames <- sampleNames[which(design==-1)]
  negGcBiases <- gcBiases[which(design==-1)]
  
  simSamplesNeg <- lapply(1:length(negSamples), function(i){
    
    simData <- varyAtacSignal(bamPath=negSamples[i], 
                              bedPath=bedPath, 
                              biasFileDir=negGcBiases[i],
                              sampleName=negSampleNames[i],
                              logFCs=logFCs,
                              effectStrength=effectStrength,
                              nFragTypes=nFragTypes,
                              fracSub=paramsGroup2$fracSub[1],
                              prob=c(paramsGroup2$prob_nf[1], 
                                     paramsGroup2$prob_mono[1], 
                                     paramsGroup2$prob_di[1], 
                                     paramsGroup2$prob_tri[1]),
                              estimateProb=paramsGroup2$estimateProb,
                              which=which,
                              minOverlap=minOverlap,
                              minGC=minGC,
                              maxGC=maxGC,
                              annotationStyle=annotationStyle,
                              genome=genome)
    
    simData$sample <- negSampleNames[i]
    simData$group <- paramsGroup2$name[1]
    simData <- as.data.table(simData)
    
    return(simData)
  })
  
  simSamplesNeg <- rbindlist(simSamplesNeg)
  simSamples <- rbind(simSamplesPos, simSamplesNeg)
  
  # even sequencing depth => make this an arguments!
  # frags dont need to be reimported over and over again
  # simSamples
  minDepth <- min(simSamples[,.N,by=sample]$N)
  simSamples <- simSamples[,.SD[sample(.N, minDepth)],by = sample]
  
  return(simSamples)
}