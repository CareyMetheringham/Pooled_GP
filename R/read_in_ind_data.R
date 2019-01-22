#' Read in the fix table
#' @param wd working directory containing fix.table
#'
#' @return an 8 col data frame:
#' SNP CROM POS ID REF ALT QUAL FILTER
#' @export
#'
#' @examples
#' read_fix_table("./extdata/test.fix")
read_fix_table <- function(fix){
  fix_table <- fread(fix)
  snp_positions <- paste(fix_table$CHROM, fix_table$POS, sep = "_")
  fix_table2 <- data.frame(snp_positions, fix_table)
  colnames(fix_table2) <-
    c("SNP", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  return(fix_table2)
}

#' Read in the genotype table
#' @param gt file working directory containing genotype.table
#'#Need to be able to vary suffix
#' @return a data frame with named rows
#' @export
#'
#' @examples
#' read_gt_table("./extdata/test.gt")
read_gt_table <- function(gt_file, suffix = "MEMq20.sorted.bam"){
  gt <- read.table(gt_file, sep = " ")
  colnames(gt) <- gsub(suffix, "", colnames(gt))
  return(gt)
}

#read in ind _info

#' Read in Individual Infomation
#'
#' @param ind_info_file path to file
#'
#' @return dataframe with two columns ID Group
#' @export
#'
#' @examples
#' read_ind_info("./extdata/example_ind_info.csv")
#' read_ind_info("./extdata/test.ind_info")
read_ind_info <- function(ind_info_file){
  info <- fread(ind_info_file)
  colnames(info) <- c("ID", "Group")
  info$ID <- paste("Individual", info$ID, sep = "")
  return(info)
}

#read in vcf data

read_vcf_file <- function(vcf_file, suffix = "MEMq20.sorted.bam"){
  ind_vcf <- read.vcfR(file=vcf_file, limit = 1e+07, cols = NULL,
                       convertNA = TRUE, checkFile = TRUE, check_keys = TRUE, verbose = TRUE)
  ind_fix <- as.data.table(getFIX(ind_vcf))

  colnames(ind_fix) <-
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  ind_fix <- cbind(paste(ind_fix$CHROM, ind_fix$POS, sep="_"), ind_fix)
  colnames(ind_fix)[1] <- "SNP"
  write.table(ind_fix,"fix.table",sep="\t",quote = FALSE, row.names = FALSE)
  ind_gt <- extract.gt(ind_vcf, element = "GT", mask = FALSE, as.numeric = TRUE,
                       return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
                       convertNA = FALSE)
  get_na <- extract.gt(ind_vcf, element = "GT", mask = FALSE, as.numeric = FALSE,
                       return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
                       convertNA = FALSE)

  #replace 0 with NA where relivant
  ind_gt[get_na=="./."] <- NA
  #remove extra stuff from names
  colnames(ind_gt) <- gsub(suffix, "", colnames(ind_gt))
  write.table(ind_gt,"genotype.table",sep="\t",quote = FALSE)
  return(list(gt = ind_gt, fix = ind_fix))
}
