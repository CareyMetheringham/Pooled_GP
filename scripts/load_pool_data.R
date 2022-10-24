#Load in data from pool_rc

#Libraries
library("argparse")
library("data.table")

######## INPUT ARGS ###########
#Define input arguments and parse
parser <-
  ArgumentParser(description = 'This reads in data from pool_rc')
parser$add_argument('--pool_rc', '-p', help = 'pool_rc file')
parser$add_argument('--info', '-i', help = 'pool info')
parser$add_argument('--snps', '-s', help = 'list of snps')
parser$add_argument('--snp_num', '-n', help = 'number of snps')
parser$add_argument('--out', '-o', help = 'Output')

xargs <- parser$parse_args()
print("args read")

######## FUNCTIONS ###########

#' Create Lists of Colnames for pools_rc
get_pool_colnames <- function(pool_info){
  colnames(pool_info) <- c("Sample", "Group", "Group2")
  mia_names <- paste("mia", pool_info$Sample, sep = "_")
  maa_names <- paste("maa", pool_info$Sample, sep = "_")
  pool_col_names <-
    c("contig",
      "pos",
      "snp",
      "P",
      "rc",
      "allele_count",
      "allele_states",
      "deletion_sum",
      "snp_type",
      "major_alleles.maa.",
      "minor_alleles.mia.",
      maa_names,
      mia_names)
  return(list(
    all = pool_col_names,
    maa = maa_names,
    mia = mia_names
  ))
}

#' Pull Top GWAS Hits from pools.rc files
find_top_snps <- function(pools_rc, hits, pool_info) {
  #get the column names
  pool_col_names <- get_pool_colnames(pool_info)$all
    colnames(pools_rc) <- pool_col_names
    pools_rc$pos <- as.character(pools_rc$pos)
    top_snps <- merge(hits, pools_rc, by = c("contig", "pos"))
  return(top_snps)
}

#' Get SNPs of Interest from GWAS Results
get_hits_from_file <- function(p_val_file, num_hits){
  hits <- fread(file = p_val_file)
  colnames(hits) <- c("contig", "pos", "snp", "P")
  hits$pos <- as.character(hits$pos)
  hits <- hits[order(as.numeric(paste(hits$P))), ]
  top_hits <- head(hits, num_hits)
  return(top_hits)
}

#' Get the Major/Minor Allele
get_allele <- function(alleles, which){
  if (which == "major") {
    pick <- 1
  }
  if (which == "minor") {
    pick <- 2
  }
  allele <-
    sapply(strsplit(as.character(alleles), "/"), function(alleles) {
      alleles <- alleles[pick]
    })
  return(allele)
}

#' Get the Major/Minor Allele Frequency
get_allele_freq <- function(pool_rc, pool_info, which){
  pool_col_names <- get_pool_colnames(pool_info)
  if (which == "major") {
    freq_data <- pool_rc[, pool_col_names$maa, with = FALSE]
  }
  if (which == "minor") {
    freq_data <- pool_rc[, pool_col_names$mia, with = FALSE]
  }
  return(as.data.frame(freq_data))
}

#' Transform allele frequencies from fraction to decimal
fraction_to_decimal <- function(frac_data) {
  dec_data <- data.frame()
  for (i in 1:length(frac_data)) {
    temp <-
      sapply(strsplit(as.character(frac_data[, i]), "/"), function(frac_data) {
        frac_data <- as.numeric(frac_data)
        frac_data[1] / frac_data[2]
      })
    dec_data <- rbind(dec_data, temp)
  }
  dec_data <- t(dec_data)
  colnames(dec_data) <- colnames(frac_data)
  return(dec_data)
}

#' Get Joined SNP ID from pools_rc format
get_snp_id <- function(pool_rc){
  snp_names <- paste(pool_rc$contig, pool_rc$pos, sep = "_")
  return(snp_names)
}

#' Exclude one named Group2 grouping
#' #can only exclude one at a time
#need to include tests for equal length
#' @param info data table with 3 columns
#' @param pool_data data output by read_in_pools_rc
#' @param group_name a group in column 3 (Group2) of the info table
#' @return pool_data excluding the named group
exclude_group2 <- function(info, pool_data, group_name){
  new_y = pool_data$y[info[, 3] != group_name]
  new_prov = pool_data$prov[info[, 3] != group_name]
  new_maa = pool_data$maa[, info[, 3] != group_name]
  new_mia = pool_data$mia[, info[, 3] != group_name]
  return(list(
    y = new_y,
    prov = new_prov,
    maa = new_maa,
    mia = new_mia,
    snp_id = pool_data$snp_id,
    major = pool_data$major,
    minor = pool_data$minor
  ))
}

#' Get random subset pools_rc
read_pools_rc_random <- function(pools_rc_file, info, hit_num){
  colnames(info) <- c("Sample", "Group", "Group2")
  random_hits <- get_random(fread(pools_rc_file, stringsAsFactors = FALSE), hit_num)
  snps_to_use <- find_top_snps(fread(pools_rc_file, stringsAsFactors = FALSE), random_hits[, 1:2], info) #ERROR HERE!!
  snp_names <- get_snp_id(snps_to_use)
  maa_freq <- get_allele_freq(snps_to_use, info, "major")
  mia_freq <- get_allele_freq(snps_to_use, info, "minor")
  maa_freq_d <- fraction_to_decimal(maa_freq)
  mia_freq_d <- fraction_to_decimal(mia_freq)
  y <- info$Group
  prov <- info$Group2
  major <- get_allele(snps_to_use$allele_states, "major")
  minor <- get_allele(snps_to_use$allele_states, "minor")
  return(list(
    y = y,
    prov = prov,
    maa = maa_freq_d,
    mia = mia_freq_d,
    snp_id = snp_names, # use as rownames
    major = major,
    minor = minor
  ))
}

#' Get random SNPs from GWAS Results
get_random <- function(pools_rc, num_hits){
  my_sample <- pools_rc[sample(nrow(pools_rc), num_hits),]
  return(my_sample)
}

######### READ IN POOL DATA ##########
info <- fread(xargs$info)
colnames(info) <- c("Sample", "Group", "Group2")
print(head(info))

#Read in top hits from the snps file
top_gwas_hits <- get_hits_from_file(xargs$snps, xargs$snp_num)
print(head(top_gwas_hits))

#Read in pool_rc
pool_rc <- fread(xargs$pool_rc, stringsAsFactors = FALSE)
head(pool_rc)
#Check that cols match those listed in get_pool_colnames and edit if not

#Pull snps from file 
snps_to_use <- find_top_snps(pool_rc, top_gwas_hits, info)
#name snps
snp_names <- get_snp_id(snps_to_use)
maa_freq <- get_allele_freq(snps_to_use, info, "major")
mia_freq <- get_allele_freq(snps_to_use, info, "minor")
maa_freq_d <- fraction_to_decimal(maa_freq)
mia_freq_d <- fraction_to_decimal(mia_freq)
y <- info$Group
prov <- info$Group2
major <- get_allele(snps_to_use$allele_states, "major")
minor <- get_allele(snps_to_use$allele_states, "minor")

poolData <- list(
    y = y,
    prov = prov,
    maa = maa_freq_d,
    mia = mia_freq_d,
    snp_id = snp_names, # use as rownames
    major = major,
    minor = minor
  )

# Save data object
save(poolData, file = xargs$out)