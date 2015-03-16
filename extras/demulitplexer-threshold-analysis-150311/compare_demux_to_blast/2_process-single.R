setwd("/data/home/em13383/other/aware-demultiplexer/extras/demulitplexer-threshold-analysis-150311/compare_demux_to_blast")

# Options
if (interactive()) {
  filen = "1out_blast_demux_merge_prob0.1_head10000.tsv.gz"
  outname = "2out_prob0.1_head10000_table.txt"
} else {
  args = commandArgs(trailingOnly = TRUE)
  filen = args[1]
  outname = args[2]
}


# Load data
df = read.table(gzfile(filen), sep="\t", header=T,
                na.strings=c("", "None", "NA"),
                stringsAsFactors=F)
df["next_score_diff"] = df$best_score - df$next_best_score
head(df)

# Remove rows with best score < 100
df2 = df[ df$best_score > 100 & !is.na(df$best_score), ]
# Remove rows where next score diff > 100
df3 = df2[ df2$next_score_diff < 100 | is.na(df2$next_score_diff), ]

# Split best hit so that it will match demux cell line
df3["best_hit_cellline"] = as.vector(sapply(df3$best_hit, function(x) { unlist(strsplit(x, " "))[2]}))

# Make table
misstab = table(df3$best_hit_cellline, df3$demux_cell_line)
write.table(misstab, file=outname, sep="\t", quote=F)
