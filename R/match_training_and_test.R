#' Find SNPs Occuring in Training and Test
#'
#' @param ees_table
#' @param gt_ind
#'
#' @return
#' @export
#'
#' @examples
#' match_snps_in_ind(fread("./extdata/test.ees_table"), read_gt_table("./extdata/test.gt"))
match_snps_in_ind <- function(ees_table, gt_ind){
  use_snps <- subset(ees_table, ees_table$SNP %in% rownames(gt_ind))
  return(use_snps)
}

#' Sort EES Table by MIA^2
#'
#' @param ees_table
#'
#' @return
#' @export
#'
#' @examples
#' order_by_ees(fread("./extdata/test.ees_table")) - test that top > bottom
order_by_ees <- function(ees_table){
  ordered_matches <-
    ees_table[order( (ees_table$EES.MIA) ^ 2, decreasing = TRUE), ]
  return(ordered_matches)
}

#' Find and Fix Major Allele Mismatches
#'
#' @param ees_table
#' @param gt_table
#' @param fix_table
#' @param pop_data
#'
#' @return
#' @export
#'
#' @examples
#' ind_fix <- read_fix_table("./extdata/test.fix")
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' test_data <- read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
#' fix_allele_mismatch(fread("./extdata/test.ees_table"), ind_gt, ind_fix, test_data)
fix_allele_mismatch <- function(ees_table, gt_table, fix_table, pop_data){
  ees_and_maa <- data.frame(ees_table, pop_data$major) #need to ensure order is correct
  colnames(ees_and_maa)[6] <- "MAJOR"
  match_by_snp <- merge(ees_and_maa, fix_table, by = "SNP") # <- This table is mostly NA!! - issue from fix?
  corrected_gt <- correct_gt(gt_table, match_by_snp)
  return(corrected_gt)
}

#' Correct Non Matching Alleles
#' @param gt
#' @param match_snps
#'
#' @return
#' @export
#'
#' @examples
correct_gt <- function(gt, match_snps){
  swapped_gt <- data.frame()
  for (i in 1:nrow(gt)){
    my_row <- rownames(gt)[i]
    if (my_row %in% match_snps$SNP){
      match_line <- match_snps[match_snps$SNP == my_row, ]
      line <- gt[i, ]
      if (match_line$MAJOR != match_line$REF){
        for (j in 1:ncol(gt)){
          if ( gt[i, j] == 0){
            line[j] <- 2
          }
          if ( gt[i, j] == 2){
            line[j] <- 0
          }
        }
      }
      swapped_gt <- rbind(swapped_gt, line)
    }
  else{next}
  }
  swapped_gt[is.na(swapped_gt)] <- 0
  return(swapped_gt)
}

#' Get a Subset of SNPs with Largest EES
#' @param ees_table
#' @param subset_size
#'
#' @return
#' @export
#'
#' @examples
#' get_ees_subset(fread("./extdata/test.ees_table"), 5)
get_ees_subset <- function(ees_table, subset_size){
  ordered_by_ees <- order_by_ees(ees_table)
  ees_subset <- head(ordered_by_ees, subset_size)
  return(ees_subset)
}

#' Subset GT Matrix by SNP list
#'
#' @param snp_list
#' @param gt
#'
#' @return
#' @export
#'
#' @examples
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' get_gt_subset(get_ees_subset(fread("./extdata/test.ees_table"), 5)$SNP, ind_gt)
get_gt_subset <- function(snp_list, gt){
  gt_subset <- subset(gt, row.names(gt) %in% snp_list) #<- THIS GIVES ERROR: rows empty
  return(gt_subset)
}

#' Match Test & Training, Correct & Subset
#'
#' @param ees_table
#' @param gt
#' @param fix
#' @param pool_data
#' @param subset_size
#'
#' @return
#' @export
#'
#' @examples
#' ind_gt <- read_gt_table("./extdata/test.gt")
#' ind_fix <- read_fix_table("./extdata/test.fix")
#' test_data <- read_in_pools_rc("./extdata/test.pool_rc", fread("./extdata/test.pool_info"), "./extdata/test.gwas", 10)
#' match_and_subset(fread("./extdata/test.ees_table"), ind_gt, ind_fix, test_data, 5)
match_and_subset <- function(ees_table, gt, fix, pool_data, subset_size){
  corrected_mismatch <-
    fix_allele_mismatch(ees_table, gt, fix, pool_data)  #<- THIS GIVES ERROR: rows empty - tested function works
  #print(head(corrected_mismatch))
  match_snps <- match_snps_in_ind(ees_table, corrected_mismatch)
  #print(head(match_snps))
  subset_snps <- get_ees_subset(match_snps, subset_size)
  #print(head(subset_snps))
  gt_subset <- get_gt_subset(subset_snps$SNP, corrected_mismatch)
  return(list(gt = gt_subset,
              ees = subset_snps))
}
