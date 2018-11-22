#' Transform allele frequencies from fraction to decimal
#' @param fracData# a data frame containing only the fractions
#'
#' @return
#' @export
#'
#' @examples
fraction_to_decimal <- function(fracData) {
  decData <- data.frame()
  for (i in 1:length(fracData)) {
    temp <-
      sapply(strsplit(as.character(fracData[, i]), "/"), function(fracData) {
        fracData <- as.numeric(fracData)
        fracData[1] / fracData[2]
      })
    decData <- rbind(temp, decData)
  }
  decData <- t(decData)
  colnames(decData) <- colnames(fracData)
  return(decData)
}
