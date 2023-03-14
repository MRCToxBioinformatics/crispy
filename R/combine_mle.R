#!/usr/bin/env Rscript
# MIT License

# Copyright (c) 2022 MRCToxBioinformatics

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(poolr))
suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
option_list <- list( 
    make_option("--significance-threshold", type='double', default=0.05, dest='sig_threshold',
        help="Threshold for significant hit [default]"),
    make_option("--ctrl-results", type='character', dest='ctrl_results',
        help="Filepath for control results"),
    make_option("--hs-results", type='character', dest='hs_results',
        help="Filepath for heat shock results"),
    make_option("--outfile", type='character',
        help="Filepath for combined results"))

opt <- parse_args(OptionParser(option_list=option_list))

add_combined_p_values <- function(obj){
  obj$stouffer <- apply(obj[,c('high.p.value', 'low.p.value')], MARGIN=1, FUN=function(x) stouffer(x)$p)

  obj$fisher <- apply(obj[,c('high.p.value', 'low.p.value')], MARGIN=1, FUN=function(x) fisher(x)$p )


  obj <- obj %>% mutate(stouffer_fdr=p.adjust(stouffer, method='BH'),
                 fisher_fdr=p.adjust(fisher, method='BH'))
  
  return(obj)
}

combine_beta <- function(obj){
  obj %>% mutate(combined_beta=(high.beta-low.beta)/2)
}

add_beta_sign_check <- function(obj){
  obj %>% mutate(beta_sign_diff=sign(high.beta)!=sign(low.beta))
}


finalise_results <- function(obj){
  obj %>%
    add_combined_p_values() %>%
    combine_beta() %>%
    add_beta_sign_check()
}


add_sig <- function(obj, sig_threshold){
  obj %>%
    mutate(sig=interaction((fisher_fdr.ctrl<sig_threshold & beta_sign_diff.ctrl),
                           (fisher_fdr.hs<sig_threshold & beta_sign_diff.hs))) %>%
    mutate(sig=recode(sig,
                      'FALSE.FALSE'='Not sig.',
                      'TRUE.FALSE'='Ctrl sig.',
                      'FALSE.TRUE'='HS sig.',
                      'TRUE.TRUE'='Both sig.')) %>%
    arrange(sig)
}

read_finalise_combine <- function(ctrl_filepath, hs_filepath, sig_threshold){
  ctrl_res <- read.table(ctrl_filepath, header=1) %>% finalise_results()
  hs_res <- read.table(hs_filepath, header=1) %>% finalise_results()
  
  keep_cols <- c('Gene', 'sgRNA', 'combined_beta', 'stouffer', 'stouffer_fdr',
                 'fisher', 'fisher_fdr', 'beta_sign_diff')
  
  combined_res <- merge(ctrl_res[,keep_cols],
                        hs_res[,keep_cols],
                        by='Gene',
                        suffixes=c('.ctrl', '.hs')) %>%
    add_sig(sig_threshold)
  
  
  return(combined_res)
}

sig_threshold <- 0.05

combined_res <- read_finalise_combine(
	opt$ctrl_results, opt$hs_results, opt$sig_threshold)

write.table(combined_res, opt$outfile, sep='\t', quote=FALSE, row.names=FALSE)

