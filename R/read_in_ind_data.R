#read in vcf data

#read in genotype table

#read in fix.table

#read in


# # 12. Read in Individual Data
# ind <- fread(ind_info)
# colnames(ind) <- c("Individual", "Health")
# ind$Individual <- paste("Individual", ind$Individual, sep = "")
# head(ind)
#
# # 13. Read in fix and genotype tables
# fix_table <- get_fix_table(wd)
# gt_table <- get_gt_table(wd)
# ind_sample_names <- gsub(".sorted.bam", "", colnames(gt_table))
# colnames(gt_table) <- ind_sample_names
