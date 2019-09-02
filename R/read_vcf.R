#read in selected rows - by gwas hit - from vcf file

#
read_vcf_file <- function(vcf_file, suffix = ""){
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
