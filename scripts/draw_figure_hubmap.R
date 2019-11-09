library(ggplot2)
library(cowplot)
library(assertthat)

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)

# Meta file for SRA ID
codedir_split <- unlist( strsplit( sub("--file=", "", args[grep("--file=", args)]), "/") )
if (length(codedir_split) > 2) {
	codedir <- paste0(codedir_split[1:(length(codedir_split)-2)], sep="/")
} else if (length(codedir_split) == 2) {
	codedir <- "./"
} else if (length(codedir_split) == 1) {
	codedir <- "../"
}
metafile <- paste0(codedir,"/data/Metadata_hubmap.txt",sep="")

# output directory
idx <- grep("--args", args)
folder <- paste0(args[idx+1], "/HumanBodyMap/", sep="")
print(paste0("Meta data file = ", metafile, sep=""))
print(paste0("Output folder = ", folder, sep=""))

# read Meta file
meta <- read.table(metafile, header=F)
colnames(meta) <- c("ID")

# read lower and upper bounds for range of optima
df <- NULL
for (ID in meta$ID) {
	print(ID)
	t_exp <- read.table(paste0(folder, ID, "/quant.sf", sep=""), header=T, sep="\t")
	t_exp$RawExp <- t_exp$NumReads / t_exp$EffectiveLength
	t_lp <- read.table(paste0(folder, ID, "/prefixgraph/salmon_lp_bound.txt", sep=""), header=F, sep="\t")
	colnames(t_lp) <- c("Name", "lb", "ub")
	t_gs <- read.table(paste0(folder, ID, "/prefixgraph/gs_maxflow_bound.txt", sep=""), header=F, sep="\t")
	colnames(t_gs) <- c("Name", "lb", "ub", "sum_gene_flow")
	# adjust the lb and ub for numerical error
	t_lp[t_lp$lb < 1e-8, "lb"] <- 0
	t_lp[t_lp$lb > t_lp$ub, "ub"] <- t_lp[t_lp$lb > t_lp$ub, "lb"]
	t_lp[abs(t_lp$lb - t_lp$ub) < 1e-8, "ub"] <- t_lp[abs(t_lp$lb - t_lp$ub) < 1e-8, "lb"]
	t_gs[t_gs$lb < 1e-8, "lb"] <- 0
	t_gs[t_gs$lb > t_gs$ub, "ub"] <- t_gs[t_gs$lb > t_gs$ub, "lb"]
	t_gs[abs(t_gs$lb - t_gs$ub) < 1e-8, "ub"] <- t_gs[abs(t_gs$lb - t_gs$ub) < 1e-8, "lb"]
	# add salmon point estimation to LP bounds and graph salmon bounds
	t_lp$point <- t_exp[match(t_lp$Name, t_exp$Name), "RawExp"]
	t_lp <- t_lp[!is.na(t_lp$point), ]
	t_gs$point <- t_exp[match(t_gs$Name, t_exp$Name), "RawExp"]
	t_gs <- t_gs[!is.na(t_gs$point), ]
	# select the one where lb != ub and salmon point estimation is within lb and ub
	tmp1 <- t_lp[t_lp$lb != t_lp$ub & t_lp$point >= t_lp$lb & t_lp$point <= t_lp$ub, c("Name", "lb", "ub", "point")]
	tmp1$Assumption <- "complete reference"
	tmp1$ID <- ID
	tmp2 <- t_gs[t_gs$lb != t_gs$ub & t_gs$point >= t_gs$lb & t_gs$point <= t_gs$ub, c("Name", "lb", "ub", "point")]
	tmp2$Assumption <- "incomplete reference"
	tmp2$ID <- ID
	if (is.null(df)) {
		df <- rbind(tmp1, tmp2)
	} else {
		df <- rbind(df, rbind(tmp1, tmp2))
	}
}

df$quantile <- (df$point - df$lb) / (df$ub - df$lb)
assert_that(all(df$ub - df$lb > 0))

# sample the incomplete reference assumption to have the same number of records as complete reference assumption
sample <- rep(TRUE, nrow(df))
ratio <- sum(df$Assumption == "complete reference") / sum(df$Assumption == "incomplete reference")
sample[df$Assumption == "incomplete reference"] <- (runif(sum(df$Assumption == "incomplete reference"), min = 0, max = 1) < ratio)

p1 <- ggplot(df[sample, ]) + geom_histogram(aes(x=quantile, fill=Assumption), alpha=0.3, position="identity") + background_grid() + scale_x_continuous(breaks=seq(0,1,0.5), labels=c("lower bound", "50%", "upper bound")) + 
	labs(title="quantile of Salmon point estimates \nwithin uncertainty range")

df2 <- df[df$point > 0, ]
df2$ratio <- df2$ub / df2$point
sample <- rep(TRUE, nrow(df2))
ratio <- sum(df2$Assumption == "complete reference") / sum(df2$Assumption == "incomplete reference")
sample[df2$Assumption == "incomplete reference"] <- (runif(sum(df2$Assumption == "incomplete reference"), min = 0, max = 1) < ratio)
print(nrow(df2[df2$Assumption == "complete reference" & abs(df2$point-df2$lb) <= 0.01*(df2$ub - df2$lb),]))
print(nrow(df2[df2$Assumption == "complete reference" & abs(df2$point-df2$ub) <= 0.01*(df2$ub - df2$lb),]))
print(nrow(df2[df2$Assumption == "complete reference" & abs(df2$point-df2$ub) > 0.01*(df2$ub - df2$lb) & abs(df2$point-df2$lb) > 0.01*(df2$ub - df2$lb),]))
p2 <- ggplot(df2[sample & df2$ratio < 10000, ]) + geom_histogram(aes(x=ratio, fill=Assumption), alpha=0.3, position="identity", bins=60) + background_grid() + 
	scale_x_log10(breaks=c(1, 2.5, 10, 100, 1000, 10000)) + labs(title="ratio between upper bound and \nSalmon point estimates")

ptmp <- plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow=1, labels="AUTO")
p <- plot_grid(ptmp, get_legend(p1), nrow=1, rel_widths=c(2,0.4))
save_plot(paste0(folder, "DeviationfromBound.pdf", sep=""), p, base_aspect_ratio = 2.4, base_height=5)

