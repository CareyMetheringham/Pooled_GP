




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
