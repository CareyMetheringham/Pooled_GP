#' Transform allele frequencies from fraction to decimal
#' @param frac_data# a data frame containing only the fractions
#'
#' @return
#' @export
#'
#' @examples
fraction_to_decimal <- function(frac_data) {
  dec_data <- data.frame()
  for (i in 1:length(frac_data)) {
    temp <-
      sapply(strsplit(as.character(frac_data[, i]), "/"), function(frac_data) {
        frac_data <- as.numeric(frac_data)
        frac_data[1] / frac_data[2]
      })
    dec_data <- rbind(temp, dec_data)
  }
  dec_data <- t(dec_data)
  colnames(dec_data) <- colnames(frac_data)
  return(dec_data)
}
