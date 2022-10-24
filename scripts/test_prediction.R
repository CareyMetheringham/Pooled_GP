#Genomic prediction from data object

#Libraries
library("argparse")
library("data.table")

######## INPUT ARGS ###########
#Define input arguments and parse
parser <-
  ArgumentParser(description = 'This tests prediction accuracy')
parser$add_argument('--es', '-e', help = 'Table of estimated effect sizes')
parser$add_argument('--info', '-i', help = 'Table of individual info')
parser$add_argument('--gt', '-g', help = 'GT table')
parser$add_argument('--out', '-o', help = 'Output file')
xargs <- parser$parse_args()

######## FUNCTIONS ###########

#' Read in the genotype table
#' @param gt file working directory containing genotype.table
#'#Need to be able to vary suffix
#' @return a data frame with named rows
#' @export
#'
#' @examples
#' read_gt_table("./extdata/test.gt")
read_gt_table <- function(gt_file, suffix = ".sorted.bam"){
  gt <- read.table(gt_file, sep = " ")
  colnames(gt) <- gsub(suffix, "", colnames(gt))
  return(gt)
}

#' Find SNPs Occuring in Training and Test
#' match_snps_in_ind(fread("./extdata/test.ees_table"), read_gt_table("./extdata/test.gt"))
match_snps_in_ind <- function(ees_table, gt_ind){
  use_snps <- subset(ees_table, ees_table$SNP %in% rownames(gt_ind))
  return(use_snps)
}

#' Subset GT Matrix by SNP list
#' @examples
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' get_gt_subset(get_ees_subset(fread("./extdata/test.ees_table"), 5)$SNP, ind_gt)
get_gt_subset <- function(snp_list, gt){
  gt_subset <- subset(gt, row.names(gt) %in% snp_list)
  return(gt_subset)
}


#' Calculate EBV from GT and ES
#' @param ees data table containing results of rrblup
#' @param gt_matrix matrix of 0, 1, 2 and NA
get_ebv <- function(ees, gt_matrix){
  gt_matrix <- gt_matrix[order(rownames(gt_matrix)), ]
  ees <- ees[order(ees$SNP), ]
  alt_matrix <- 2 - gt_matrix
  effect_mia <-  t(gt_matrix) %*% ees$EES.MIA
  effect_maa <-  t(alt_matrix) %*% ees$EES.MAA
  ebv <- effect_mia + effect_maa
  ebv <- as.vector(ebv)
  names(ebv) <- colnames(gt_matrix)
  return(ebv)
}

#' Create EBV Table
#' @param ind_info phenotypic assignment of test individuals
#' @param ebv estimated breeding value
#' @return table of estimated breeding values and assigned groups for test individuals
#' @examples
create_ebv_table <- function(ind_info, ebv){
  ebv_df <- data.frame(ID = names(ebv), EBV = as.numeric(ebv))
  ebv_table <- merge(ind_info, ebv_df, by = "ID")
  return(ebv_table)
}

######## TEST PREDICTION ACCURACY ###########

#Get effect size table
ees_table <- fread(xargs$es)

#Get GT matrix 0,1,2 or -1, 0, 1
#strip .bam if still present
gt <- read_gt_table(xargs$gt)

#Get subset of matched snps
matched_ees <- match_snps_in_ind(ees_table, gt)
matched_gt <- get_gt_subset(matched_ees$SNP, gt)

print(head(matched_gt))

#Estimate breeding values
ebv <- get_ebv(matched_ees, matched_gt)

#Get individual information
ind_info <- fread(xargs$info)

#Write GEBV to output file
ebv_table <- create_ebv_table(ind_info = ind_info, ebv = ebv)
write.csv(ebv_table, file = xargs$out, colnames = F)