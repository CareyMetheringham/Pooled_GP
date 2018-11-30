#' Read in the fix table
#' @param wd working directory containing fix.table
#'
#' @return an 8 col data frame:
#' SNP CROM POS ID REF ALT QUAL FILTER
#' @export
#'
#' @examples
#' read_fix_table("./extdata")
read_fix_table <- function(wd){
  fix <- paste(wd, "fix.table", sep = "/")
  fix_table <- fread(fix)
  snp_positions <- paste(fix_table$CHROM, fix_table$POS, sep = "_")
  fix_table2 <- data.frame(snp_positions, fix_table)
  colnames(fix_table2) <-
    c("SNP", "CROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  return(fix_table2)
}

#' Read in the genotype table
#' @param wd working directory containing genotype.table
#'
#' @return a data frame with named rows
#' @export
#'
#' @examples
#' read_gt_table("./extdata")
read_gt_table <- function(wd){
  genotype <- paste(wd, "genotype.table", sep = "/")
  genotype_table <- read.table(genotype, sep = "\t")
  snps_in_ind <- genotype_table$V1
  gt <- as.data.frame(genotype_table[, -1])
  colnames(gt) <- gsub(".sorted.bam", "", colnames(gt))
  rownames(gt) <- snps_in_ind
  return(gt)
}

#read in ind _info

#read in vcf data
