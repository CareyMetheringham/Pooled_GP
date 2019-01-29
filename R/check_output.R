#' Calculate Accuracy of Predictions
#'
#' @param ebv_table
#'
#' @return
#' @export
#'
#' @examples
calculate_accuracy <- function(ebv_table, percentage = 1, G1val = 1, G2val = 2){
  ebv_sort <- ebv_table[order(ebv_table$EBV), ]
  group_size <- nrow(ebv_table) / 2 * percentage
  G1 <- head(ebv_sort, group_size)
  G2 <- tail(ebv_sort, group_size)
  G1count <- count(G1$Group)
  G2count <- count(G2$Group)
  G1correct <- G1count$freq[G1count$x == G1val]
  G2correct <- G2count$freq[G2count$x == G2val]
  accuracy <- (G1correct + G2correct) / (group_size * 2)
  return(accuracy)
}

#' Title
#'
#' @param ind_info
#' @param ebv
#'
#' @return
#' @export
#'
#' @examples
create_ebv_table <- function(ind_info, ebv){
  ebv_df <- data.frame(ID = names(ebv), EBV = as.numeric(ebv))
  ebv_table <- cbind(ind_info, ebv_df)
  return(ebv_table)
}

#' Title
#'
#' @param ebv_table
#'
#' @return
#' @export
#'
#' @examples
calculate_correlation <- function(ebv_table){
  correlation <-
    cor(ebv_table$Group, ebv_table$EBV, method = "pearson")
  return(correlation)
}
