library("ggplot2")

setwd("/data/home/em13383/other/aware-demultiplexer/extras/demulitplexer-threshold-analysis-150311/blast-multiplexed-seqs")

# Load data
data = read.table(gzfile("4out_blast-output-parsed.tsv.gz"), sep="\t", header=T,
                  na.strings="None")
head(data)

if (!isinteractive()) {
  png()
}

# If there is a next best score then calc score difference
data["next_best_diff"] = data$best_score - data$next_best_score

# Make plot of best score by cell line
ggplot(data, aes(x=best_score, fill=best_hit)) +
  geom_histogram(position="identity", alpha=0.5)
# Make plot of next best score diff by cell line
ggplot(data, aes(x=next_best_diff, fill=best_hit)) +
  geom_histogram(position="identity", alpha=0.5)

# Print how many hits there are for each cell line
total.num = nrow(data)
num.hits = data.frame(table(data$best_hit))
num.hits["Proportion"] = num.hits$Freq / total.num
num.hits
total.num


