arg <- commandArgs(trailingOnly = TRUE)
# arg[1] = reference BED file
# arg[2] = output file
# arg[3..N] = input .tab files
ref_bed <- arg[1]
out_file <- arg[2]
tab_ls <- arg[3:length(arg)]

# Extract locus_tag from reference BED (column 4)
locus_tag <- read.delim(ref_bed, header = FALSE, sep = "\t")$V4

# Derive sample names from filenames
tab_nm <- gsub("\\..*\\.tab$", "", basename(tab_ls))

# Read all .tab files and rename value column to sample name
tab_list <- mapply(function(x, y) {
  df <- read.delim(x, header = TRUE, sep = "\t", check.names = FALSE)
  colnames(df)[2] <- y
  return(df)
}, tab_ls, tab_nm, SIMPLIFY = FALSE)
tab_list <- lapply(tab_list, function(x) subset(x, select = -c(locus_tag)))
tab_df <- do.call(cbind, tab_list)
tab_df$locus_tag <- locus_tag
tab_df <- tab_df[, c("locus_tag", head(colnames(tab_df), -1))] # Reorder columns
write.table(tab_df, out_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
