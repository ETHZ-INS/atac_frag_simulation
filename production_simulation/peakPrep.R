library(chromVAR)
library(memes)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(MotifDb)
library(motifmatchr)
library(Matrix)
library(TFBSTools)
library(universalmotif)

getNonRedundantMotifs <- function(format=c("PFMatrix","universal","PWMatrix"),
                                  species=c("Hsapiens","Mmusculus")){
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb, c(species,"HOCOMOCO"))
  pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
  modf <- data.frame(row.names=names(motifs),
                     TF=gsub(pat,"",names(motifs)),
                     grade=gsub(".+\\.","",names(motifs)))
  modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
  modf <- modf[!duplicated(modf$TF),]
  motifs <- motifs[row.names(modf)]
  switch(match.arg(format),
         universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
         PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
           universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
           modf$TF)),
         PWMatrix=do.call(TFBSTools::PWMatrixList, 
                          setNames(universalmotif::convert_motifs(motifs, 
                                                                  class="TFBSTools-PWMatrix"), modf$TF))
  )
}

getpmoi <- function(genome,
                    peakpath,
                    outDir,
                    spec=c("Hsapiens","Mmusculus"),
                    seqStyle=c("ensembl","UCSC")){
  # choose the file that contains the correct names from HOCOMOCO v11
  if (spec=="Hsapiens") {
    motifnames <- fread("/mnt/plger/fgerbaldo/BenchmarkTFactivity/BMScripts/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv")
  } else if (spec=="Mmusculus") {
    motifnames <- fread("/mnt/plger/fgerbaldo/BenchmarkTFactivity/BMScripts/HOCOMOCOv11_core_annotation_MOUSE_mono.tsv")
  }
  # Read in and resize the peaks
  peakfile <- file(peakpath)
  peaks <- chromVAR::getPeaks(peakfile, 
                              sort_peaks = TRUE)
  peaks <- resize(peaks, 
                  width = 300, 
                  fix = "center")
  peaks <- keepStandardChromosomes(peaks, 
                                   pruning.mode = "coarse")
  
  seqlevelsStyle(genome) <- seqStyle
  peak_seqs <- get_sequence(peaks, genome)
  
  # Get the motifs in universal format required by memes
  motifs <- getNonRedundantMotifs("universal", species = spec)
  
  BANP_motif <- readRDS("/mnt/plger/datasets/Grand2021_BANP_ATAC_GSE155601/BANP.PFMatrix.rds")
  BANP_universalmotif <- convert_motifs(BANP_motif, class = "universalmotif-universalmotif")
  motifs$BANP <- BANP_universalmotif
  
  # Correct inconvential names for motifs
  for (i in seq_along(motifs)){ 
    for (j in seq_along(motifnames$Model)) {
      if (motifs[[i]]@altname == motifnames$Model[[j]])
        names(motifs)[[i]] <- motifnames$`Transcription factor`[[j]]
    }
  }
  
  # Obtain the positions of motif instances which are later required as input for runATAC
  pmoi <- runFimo(peak_seqs, 
                  motifs, 
                  meme_path="/common/meme/bin/", 
                  skip_matched_sequence=TRUE)
  saveRDS(pmoi, file.path(outDir, "pmoi.rds"))
}


# Peak filtering
peaks <- fread("/mnt/plger/datasets/ENCODE_lymphoblastoid_ATAC/peaks/mergedSamples_atac_peaks.bed_peaks.narrowPeak", col.names=c("chr", "start", "end"), select=1:3)
peaks <- makeGRangesFromDataFrame(peaks)
bl <- rtracklayer::import("/mnt/reference/reference/Homo_sapiens/hg38.blacklist.chrM.bed")
seqlevelsStyle(bl) <- "UCSC"

peaks_filtered <- subsetByOverlaps(peaks, bl, invert=TRUE)
rtracklayer::export.bed(peaks_filtered, con='./atac_peaks/filtered_mergedSamples_atac_peaks.bed_peaks.narrowPeak')

# Finding motifs within peaks
genome <- BSgenome.Hsapiens.UCSC.hg38
peakpath <- './atac_peaks/filtered_mergedSamples_atac_peaks.bed_peaks.narrowPeak'
getpmoi(genome, peakpath, outDir="./atac_peaks",spec = "Hsapiens",seqStyle = "UCSC")
