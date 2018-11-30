match_snps_in_ind <- function(ees_table, gt_ind){
  use_snps <- subset(ees_table, ees_table$SNP %in% row.names(gt_ind))

}

order_by_ees <- function(ees_table){
  ordered_matches <- ees_table[order((ees_table$EES.MIA) ^ 2, decreasing = TRUE),]
}

fix_allele_mismatch <- function(ees_table, gt_table, fix_table, pop_data){
  ees_and_maa <- data.frame(ees_table, pop_data$major)
  colnames(ees_and_maa)[6] <- "MAJOR"
  match_by_snp <- merge(ees_and_maa, fix_table, by = "SNP")
  corrected_gt <- correct_gt(gt_table, match_by_snp)
  return(corrected_gt)
}

#' Correct Non Matching alleles
#' @param gt
#' @param match_snps
#'
#' @return
#' @export
#'
#' @examples
correct_gt <- function(gt, match_snps){
  swapped_gt <- gt
  for (i in 1:nrow(gt)){
    my_row <- rownames(gt)[i]
    if (my_row %in% row.names(match_snps)){
      if (match_snps[my_row, "MAJOR"] != match_snps[my_row, "REF"]){
        # Mismatch - swap values
        swapped_gt[i, ] <- sub("0", "NO", gt[i, ])
        swapped_gt[i, ] <- sub("2", "0", gt[i, ])
        swapped_gt[i, ] <- sub("NO", "2", gt[i, ])
      }
    }
  }
  swapped_gt[is.na(swapped_gt)] <- 0
  return(swapped_gt)
}

get_ees_subset <- function(ees_table, subset_size){
  ordered_by_ees <- order_by_ees(ees_table)
  ees_subset <- head(ordered_by_ees, subset_size)
  return(ees_subset)
}

get_gt_subset <- function(snp_list, gt){
  gt_subset <- subset(gt, row.names(gt) %in% snp_list)
}

match_and_subset <- function(ees_table, gt, fix, pool_data, subset_size){
  corrected_mismatch <- fix_allele_mismatch(ees_table, ind_gt, ind_fix, pool_data)
  match_snps <- match_snps_in_ind(ees_table, corrected_mismatch)
  subset_snps <- get_ees_subset(match_snps, subset_size)
  gt_subset <- get_gt_subset(subset_snps$SNP, corrected_mismatch)
  return(list(gt = gt_subset,
              ees = subset_snps))
}

#' Find the Subset of Genotypes
#' #' @param gt
#' #' @param snps
#' #' @param match_table
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' get_gt_subset <- function(gt, snps, match_table){
#'   top_snps <- head(match_table, as.numeric(snps))
#'   order_top_snps <- top_snps[order(rownames(top_snps)), ]
#'   rows_to_use <- top_snps$SNP
#'   gt_subset <- subset(gt, row.names(gt) %in% rows_to_use)
#'   order_gt_subset <- gt_subset[order(rownames(gt_subset)), ]
#'   return(list(gt = order_gt_subset,
#'               data = order_top_snps))
#}


#need to get and store major and minor alleles


# # 14. Find the matching SNPs
# use_ind <- subset(ees_table, row.names(ees_table) %in% row.names(gt_table))
# match_snps <-
#   merge(as.data.frame(fix_table), as.data.frame(use_ind), by = "SNP")
#T
# # 15. Sort the data table by EES
# sorted_by_ees <-
#   match_snps[order((match_snps$MIA.EES) ^ 2, decreasing = TRUE),]
# print(head(sorted_by_ees, 10))
# #count total number
# total_snps_matched <- nrow(sorted_by_ees)
#
# # 16. Correct none matching SNPs
# corrected_gt <- correct_non_matching(gt_table, sorted_by_ees)
