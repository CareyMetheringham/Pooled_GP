matched_ees <- match_snps_in_ind(ees_table, test_data$gt)
matched_gt <- get_gt_subset(matched_ees$SNP, test_data$gt)


print("Estimate breeding values")
ebv <- get_ebv(matched_ees, matched_gt)
print("Correlation of EBV and true BV in test population")
print(cor(ebv, as.vector(test_data$bv)))
plot(as.vector(test_data$bv),
     ebv,
     xlab = "Input Breeding Value",
     ylab = "Estimated Breeding Value")
print("Correlation of EBV and observed phenotype in test population")
print(cor(ebv, as.vector(test_data$ph)))
plot(as.vector(test_data$ph),
     ebv,
     xlab = "Phenotypic Value",
     ylab = "Estimated Breeding Value")
}