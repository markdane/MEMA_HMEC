#!/usr/bin/env Rscript

#author: "Mark Dane"
# 1/15/19

library(MEMA)#merge, annotate and normalize functions
library(parallel)#use multiple cores for faster processing
library(RUVnormalize)
library(ruv)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
library(rrscale)

######
rrscale_vector <- function(x, ...) {
  x_rr <- rrscale(as.matrix(x), ...)
  return(x_rr$RR)
}
######

studyName <- "hmec240l_ss4"
inputPath <- "/Users/dane/Documents/MEMATechPaper/HMEC240L_HGF/Data/HMECData/"
ofname <-  paste0("/Users/dane/Documents/MEMATechPaper/HMEC240L_HGF/Data/HMECData/",studyName,"_rr_RUV_level3.tsv")
k <- 256
verbose <- TRUE

#TODO
#Load in level3 data
#Filter to columns of interest
#Untransform logged data to create l2 dataset


#Read the annotated data for all plates in the study
l2 <- dir(inputPath,
          recursive = TRUE,
          pattern = "Level2",
          full.names = TRUE) %>%
  map(read_tsv) %>%
  bind_rows()

l2_rr <- l2 %>%
  select(Spot_PA_SpotCellCount,
         Cytoplasm_PA_Gated_KRT19PositiveProportion,
         Cytoplasm_PA_Gated_KRT19Positive_SpotCellCount,
         Nuclei_PA_Gated_EdUPositiveProportion,
         Cells_PA_Gated_EdUKRT19PositiveProportion) %>%
  transmute_all(list(RR = rrscale_vector), zeros = .001) %>%
  bind_cols(l2) %>%
  mutate(BW = paste0(Barcode,"_",Well),
         k = k) %>%
  data.table

signalsMinimalMetadata <- grep("_SE",
                               grep("_RR|Barcode|^Well$|^Spot$|^PrintSpot$|^Ligand$|^ECMp$|^Drug$|^Drug1Conc$|^ArrayRow$|^ArrayColumn$|^CellLine$",
                                    colnames(l2_rr), value=TRUE), 
                               value=TRUE, invert=TRUE)

#RUVLoess normalize all the rr signals
  if(verbose)  message(paste("Normalizing", studyName,"\n"))
  nDT <- norm_RR_RUVResiduals(l2_rr[,signalsMinimalMetadata, with = FALSE], k) %>%
    mutate(NormMethod = "RR_RUVResiduals") %>%
    right_join(l2_rr, by = c("BW", "PrintSpot")) %>%
    mutate_if(is.numeric, signif(6))


if (verbose) message(paste("Writing level 3 file to disk\n"))
write_csv(nDT, ofname)

message(paste("Elapsed time to normalize ",studyName, Sys.time()-startTime, "\n"))


