library(ggplot2)
library(cowplot)
library(DESeq2)
library(tximport)

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
codedir_split <- strsplit(args[1], "/")
codedir <- paste0(codedir_split[1:(length(codedir_split)-2)], sep="/")
Dir <- args[2]

####################################### summarizing bounds #######################################
t <- read.table(paste0(Dir,"/GEUVADIS/salmon_Full_ERR188021/prefixgraph/gs_bound_round0.txt",sep=""), header=F, sep="\t")
colnames(t) <- c("Name", "lb", "ub", "sum_gene_flow")
t <- t[t$ub > 0, ]
t[t$lb < 0, "lb"] <- 0
t[t$ub < 0, "ub"] <- 0
t[t$ub < t$lb,  "ub"] <- t[t$ub < t$lb,  "lb"]
t$ratio_ub <- (t$ub - t$lb) / t$ub

genesmap <- read.table(paste0(Dir,"/gencode.v26.Gene_Trans_Map.txt",sep=""))
colnames(genesmap) <- c("GeneID", "TransID")
count <- data.frame(table(genesmap$GeneID))

t$GeneID <- genesmap[match(t$Name, genesmap$TransID), "GeneID"]
t$NumTrans <- count[match(t$GeneID, count$Var1), "Freq"]
p1 <- ggplot(t) + geom_histogram(aes(x=ratio_ub)) + scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) + background_grid() + labs(x="Gap ratio", y="Count")

count <- data.frame(table(t$NumTrans))
df <- data.frame(NumTrans=rep(as.numeric(count$Var1),4), GapRatio=c(rep(c("[0,0.1)","[0.1,0.5)","[0.5,0.9)","[0.9,1]"), nrow(count))), Count=rep(0, nrow(count)*4))
for (i in seq(nrow(df))) {
	if (df[i, "GapRatio"] == "[0,0.1)") {
		df[i, "Count"] <- nrow(t[t$NumTrans==df[i, "NumTrans"] & t$ratio_ub>=0 & t$ratio_ub<0.1, ])
	} else if (df[i, "GapRatio"] == "[0.1,0.5)") {
		df[i, "Count"] <- nrow(t[t$NumTrans==df[i, "NumTrans"] & t$ratio_ub>=0.1 & t$ratio_ub<0.5, ])
	} else if (df[i, "GapRatio"] == "[0.5,0.9)") {
		df[i, "Count"] <- nrow(t[t$NumTrans==df[i, "NumTrans"] & t$ratio_ub>=0.5 & t$ratio_ub<0.9, ])
	} else {
		df[i, "Count"] <- nrow(t[t$NumTrans==df[i, "NumTrans"] & t$ratio_ub>=0.9, ])
	}
}
df$Count <- df$Count / count[match(df$NumTrans, count$Var1), "Freq"]
p2 <- ggplot(df) + geom_bar(aes(x=NumTrans, y=Count, fill=GapRatio), stat="identity") + scale_x_continuous(limits=c(0,40)) + labs(x="Number annotated transcripts", y="Percentage") + scale_fill_brewer(palette="Set2")

p <- plot_grid(p1, p2, labels="AUTO", rel_widths=c(0.45, 0.55))
save_plot(paste0(Dir,"/ERR188021_abundance_bound.pdf",sep=""), p, base_aspect_ratio = 2.2)


####################################### DE detection and example #######################################
ReadMeta <- function(metafile) {
	Centers  <- c()
	IDs <- c()
	con = file(metafile, "r")
	line <- readLines(con, n = 1)
	while (length(line) > 0) {
		strs <- strsplit(line, "\t")[[1]]
		Centers <- c(Centers, strs[2])
		IDs <- c(IDs, strs[length(strs)])
		line <- readLines(con, n = 1)
	}
	close(con)
	t_meta <- data.frame(Center=Centers, ID=IDs)
	return(t_meta)
}

t_meta <- ReadMeta(paste0(codedir,"/data/Metadata.txt",sep=""))
C <- unique(t_meta$Center)

# DESeq2 to infer DE transcripts
files <- paste0(Dir,"/GEUVADIS/salmon_Full_", t_meta$ID, "/quant.sf", sep="")
names(files) <- t_meta$ID
txi <- tximport(files, type = "salmon", txOut = TRUE)
dds <- DESeqDataSetFromTximport(txi, t_meta, ~Center)
dds <- DESeq(dds)
res <- results(dds, alpha=0.01)
print(sum(res$padj < 0.01, na.rm=TRUE))
res2 <- res[!is.na(res$padj),]

genesmap <- read.table(paste0(Dir,"/gencode.v26.Gene_Trans_Map.txt",sep=""))
colnames(genesmap) <- c("GeneID", "TransID")
count <- data.frame(table(genesmap$GeneID))
res2$GeneID <- genesmap[match(rownames(res2), genesmap$TransID), "GeneID"]
res2$NumTrans <- count[match(res2$GeneID, count$Var1), "Freq"]

t_ard <- read.table(paste0(Dir,"/GEUVADIS/bounds_meanflowARD.txt",sep=""), header=F, sep="\t")
colnames(t_ard) <- c("GeneID", "MeanARD")
res2$gene_ard <- t_ard[match(res2$GeneID, t_ard$GeneID), "MeanARD"]

t_gap_overlap <- read.table(paste0(Dir, "GEUVADIS/bounds_weightedJaccard.txt",sep=""), header=F, sep="\t")
colnames(t_gap_overlap) <- c("Name", "gap_overlap")
res2$gap_overlap <- t_gap_overlap[match(rownames(res2), t_gap_overlap$Name), "gap_overlap"]
df <- data.frame(res2[res2$padj < 0.01 & res2$gene_ard < 0.2, ])
df$NumTrans <- factor(df$NumTrans, level=seq(1,5))
df <- df[!is.na(df$NumTrans),]
p1 <- ggplot(df) + geom_violin(aes(x=NumTrans, y=gap_overlap)) + background_grid() + 
	labs(y="Similarity", x="Number annotated transcripts", title="Weighted Jaccard similarity of bounds between centers")

t_dist_example <- read.table(paste0(Dir,"GEUVADIS/bounds_example_dist_ENST00000422247.6.txt",sep=""), header=F, sep="\t")
colnames(t_dist_example) <- c("Name","Center","position","count")
p22 <- ggplot(t_dist_example) + geom_step(aes(x=position, y=count, color=Center), alpha=0.5, size=1) + background_grid() + labs(x="Possible abundance within bounds", y="Count of samples", title="Density of bounds of ENST00000422247.6")

tmp_example <- read.table(paste0(Dir,"GEUVADIS/bounds_exampleENST00000375754.8.txt",sep=""), header=F, sep="\t")
colnames(tmp_example) <- c("Name", "Sample", "Center", "lb", "ub", "TPM")
t_example <- rbind(tmp_example[,c(1,2,3,6)], tmp_example[,c(1,2,3,6)])
t_example$bound <- c(tmp_example$lb, tmp_example$ub)
t_example$bound_dir <- c(rep("lower bound",nrow(tmp_example)), rep("upper bound",nrow(tmp_example)))
p3 <- ggplot(t_example) + geom_line(aes(x=bound, y=TPM, color=Center, group=Sample), alpha=0.5) + geom_point(aes(x=bound, y=TPM, color=Center, shape=bound_dir), size=2, alpha=0.7) + background_grid() +
	labs(x="Possible abundance within bounds", y="Salmon TPM", title="Bounds of ENST00000375754.8", color="sequencing center", shape="")

t_dist_example <- read.table(paste0(Dir,"GEUVADIS/bounds_example_dist_ENST00000375754.8.txt",sep=""), header=F, sep="\t")
colnames(t_dist_example) <- c("Name","Center","position","count")
p33 <- ggplot(t_dist_example) + geom_step(aes(x=position, y=count, color=Center), alpha=0.5, size=1) + background_grid() + labs(x="Possible abundance within bounds", y="Count of samples", title="Density of bounds of ENST00000375754.8")

t_dist_example <- read.table(paste0(Dir,"GEUVADIS/bounds_example_dist_ENST00000265187.4.txt",sep=""), header=F, sep="\t")
colnames(t_dist_example) <- c("Name","Center","position","count")
p44 <- ggplot(t_dist_example) + geom_step(aes(x=position, y=count, color=Center), alpha=0.5, size=1) + background_grid() + labs(x="Possible abundance within bounds", y="Count of samples", title="Density of bounds of  ENST00000265187.4")

p <- ggdraw() + 
	draw_plot(p1, 0.25, 0.68, 0.5, 0.32) + 
	draw_plot(p3 + theme(legend.position="none"), 0.01, 0.36, 0.49, 0.32) + 
	draw_plot(p33 + theme(legend.position="none"), 0.5, 0.36, 0.49, 0.32) +
	draw_plot(p22 + theme(legend.position="none"), 0.01, 0.04, 0.49, 0.32) + 
	draw_plot(p44 + theme(legend.position="none"), 0.5, 0.04, 0.49, 0.32) +
	draw_plot(get_legend(p3 + theme(legend.position="bottom")), 0.2, 0, 0.8, 0.04) +
	draw_plot_label(c("A", "B", "C", "D", "E"), c(0.24,0,0.5,0,0.5), c(1,0.68,0.68,0.36,0.36))
save_plot(paste0(Dir,"abundance_gap_example.pdf",sep=""), p, base_aspect_ratio = 1, base_height = 10)
