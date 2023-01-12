#!/usr/bin/env Rscript
# 12.12.2018
# Daniel Schreyer

#----------------------------------------------------------------
# Correlation between ref allele freq and target allele freq Plot
# Takes Imputation info file and target file as an Input and outputs two differently colored scatter plots
# one plot is colored based on the SNP r-squared values 
# and one based on the difference between the allele frequencies
#----------------------------------------------------------------

# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)

freq_comp_plot <- function(freq, plot) {

  info <- read.table(freq, header = F, sep ="\t")
  names(info) <- c('DatasetMAF', 'RefPAnelMAF', 'INFO')

  INFO.thresh = 0.8
  info <- filter(info, INFO != "-" | !is.na(INFO)) %>% filter(INFO > INFO.thresh)

  # filter every Nth SNP to extract [] SNPs [Default = 20000]
  # subset = 200000
  # info <- info[1:subset, ]

  plot.rsq.colored <- ggplot(info , aes(x = DatasetMAF, y = RefPAnelMAF, color = INFO)) +
    geom_point(size = 0.9) + theme_classic() + 
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
    labs(x = "Ref Allele Frequency (Uploaded Samples)", 
        y = "Ref Allele Frequency (Reference Panel)") +
    geom_text(color = "black", x = 0.2, y = 0.95, 
              label = paste(nrow(info),"SNPs", sep = " ")) +
    scale_color_gradient(low = "lightblue", high = "darkblue")
  ggsave(filename = plot, plot = plot.rsq.colored, width = 7, height = 7, units = "in")
  return(plot)
}

freq_comp_plot1 <- function(freq, plot) {

  info <- read.table(freq, header = F, sep ="\t")
  names(info) <- c('DatasetMAF', 'RefPAnelMAF', 'INFO')

  INFO.thresh = 0.8
  info <- filter(info, INFO != "-" | !is.na(INFO)) %>% filter(INFO > INFO.thresh)

  # filter every Nth SNP to extract [] SNPs [Default = 20000]
  subset = 200000
  info <- info[1:subset, ]

  plot.rsq.colored <- ggplot(info , aes(x = DatasetMAF, y = RefPAnelMAF, color = INFO)) +
    geom_point(size = 0.9) + theme_classic() + 
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) + 
    labs(x = "Ref Allele Frequency (Uploaded Samples)", 
        y = "Ref Allele Frequency (Reference Panel)") +
    geom_text(color = "black", x = 0.2, y = 0.95, 
              label = paste(nrow(info),"SNPs", sep = " ")) +
    scale_color_gradient(low = "lightblue", high = "darkblue")
  ggsave(filename = plot, plot = plot.rsq.colored, width = 7, height = 7, units = "in")
  return(plot)
}

freq_comp_plot1("sanger_afr.chr_all.freq.final.txt", "sanger_afr.chr_all.freq_2M.tiff")
freq_comp_plot1("sanger_kgp.chr_all.freq.final.txt", "sanger_kgp.chr_all.freq_2M.tiff")
freq_comp_plot1("sanger_hrc.chr_all.freq.final.txt", "sanger_hrc.chr_all.freq_2M.tiff")

freq_comp_plot("sanger_afr.chr_all.freq.final.txt", "sanger_afr.chr_all.freq.tiff")
freq_comp_plot("sanger_kgp.chr_all.freq.final.txt", "sanger_kgp.chr_all.freq.tiff")
freq_comp_plot("sanger_hrc.chr_all.freq.final.txt", "sanger_hrc.chr_all.freq.tiff")