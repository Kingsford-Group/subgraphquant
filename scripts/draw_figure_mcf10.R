library(ggplot2)
library(cowplot)
library(assertthat)
library(DESeq2)
library(tximport)

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
metafile <- paste0(codedir,"/data/Metadata_mcf10.txt",sep="")

# output directory
idx <- grep("--args", args)
folder <- paste0(args[idx+1], "/MCF10/", sep="")
print(paste0("Meta data file = ", metafile, sep=""))
print(paste0("Output folder = ", folder, sep=""))

# read Meta file
meta <- read.table(metafile, header=F, sep="\t")
meta <- meta[, c(5, 8, 10)]
colnames(meta) <- c("ID", "treatment", "time")

# read IValue
iv <- read.table(paste0(folder, "IValue_mean_ext.txt", sep=""), header=F, sep="\t")
colnames(iv) <- c("Name", "Group1Mean", "Group2Mean", "IV", "IV25")

# DE detection on transcript level
files <- paste0(folder, meta$ID, "/quant.sf", sep="")
names(files) <- meta$ID
txi <- tximport(files, type = "salmon", txOut = TRUE)
dds <- DESeqDataSetFromTximport(txi, meta, ~treatment)
dds <- DESeq(dds)
res <- results(dds, alpha=0.01)
print(sum(res$padj < 0.01, na.rm=TRUE))
res <- res[!is.na(res$padj),]
res$IV <- iv[match(rownames(res), iv$Name), "IV"]
res$IV25 <- iv[match(rownames(res), iv$Name), "IV25"]

# plot histogram
p2 <- ggplot(as.data.frame(res[res$padj < 0.01, ])) + geom_histogram(aes(x=1-IV25)) + background_grid() + labs(x = "reference completeness", title = "Number of DE transcripts that are unreliable \nunder the reference completeness parameter")

# draw example curve
t <- read.table(paste0(folder, "IValue_curve_example.txt", sep=""), header=F, sep="\t")
colnames(t) <- c("Name", "Group", "Sample", "reference_proportion", "lb", "ub")
df <- data.frame(Name=rep(t$Name, 2), Group=rep(t$Group, 2), Sample=rep(t$Sample, 2), reference_proportion=rep(t$reference_proportion, 2), Expression=c(t$lb, t$ub), Type=c(rep("lower bound", nrow(t)), rep("upper bound", nrow(t))))

p3 <- ggplot(df[df$Name == "ENST00000010404.6" & df$Sample == "Mean", ]) + geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), color=Group, linetype=Type), size=1.2) + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000010404.6" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000010404.6" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "Expression"], ub=df[df$Name == "ENST00000010404.6" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="EGF Stimulation + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#F8766D") + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000010404.6" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000010404.6" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "Expression"], ub=df[df$Name == "ENST00000010404.6" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="no EGF + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#00BFC4") + 
	geom_vline(xintercept=1-0.0559912296325953, color="black") + geom_vline(xintercept=1-0.109375733685513, color="red") + background_grid() + labs(title = "Ranges of optima of ENST00000010404.6", x = "reference completeness", y = "normalized abundance") +
	theme(legend.title=element_blank())
p4 <- ggplot(df[df$Name == "ENST00000396832.5" & df$Sample == "Mean", ]) + geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), color=Group, linetype=Type), size=1.2) + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000396832.5" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000396832.5" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "Expression"], ub=df[df$Name == "ENST00000396832.5" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="EGF Stimulation + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#F8766D") + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000396832.5" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000396832.5" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "Expression"], ub=df[df$Name == "ENST00000396832.5" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="no EGF + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#00BFC4") + 
	geom_vline(xintercept=1-0.0426357585576571, color="black") + geom_vline(xintercept=1-0.0641238355024918, color="red") + background_grid() + labs(title = "Ranges of optima of ENST00000396832.5", x = "reference completeness", y = "normalized abundance") +
	theme(legend.title=element_blank())
p5 <- ggplot(df[df$Name == "ENST00000377482.9" & df$Sample == "Mean", ]) + geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), color=Group, linetype=Type), size=1.2) + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000377482.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000377482.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "Expression"], ub=df[df$Name == "ENST00000377482.9" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="EGF Stimulation + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#F8766D") + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000377482.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000377482.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "Expression"], ub=df[df$Name == "ENST00000377482.9" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="no EGF + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#00BFC4") + 
	geom_vline(xintercept=1-0.655461069604338, color="black") + geom_vline(xintercept=1-0.694297745102267, color="red") + background_grid() + labs(title = "Ranges of optima of ENST00000377482.9", x = "reference completeness", y = "normalized abundance") +
	theme(legend.title=element_blank())
p6 <- ggplot(df[df$Name == "ENST00000337393.9" & df$Sample == "Mean", ]) + geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), color=Group, linetype=Type), size=1.2) + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000337393.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000337393.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="EGF Stimulation + DMSO", "Expression"], ub=df[df$Name == "ENST00000337393.9" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="EGF Stimulation + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#F8766D") + 
	geom_ribbon(data=data.frame(reference_proportion=df[df$Name == "ENST00000337393.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "reference_proportion"], lb=df[df$Name == "ENST00000337393.9" & df$Sample == "Mean" & df$Type=="lower bound" & df$Group=="no EGF + DMSO", "Expression"], ub=df[df$Name == "ENST00000337393.9" & df$Sample == "Mean" & df$Type=="upper bound" & df$Group=="no EGF + DMSO", "Expression"]), aes(x=1-reference_proportion, ymin=lb, ymax=ub), alpha=0.2, fill="#00BFC4") + 
	background_grid() + labs(title = "Ranges of optima of ENST00000337393.9", x = "reference completeness", y = "normalized abundance") +
	theme(legend.title=element_blank())

ptmp <- plot_grid(p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), p5 + theme(legend.position="none"), p6 + theme(legend.position="none"), p2, nrow=2, labels="AUTO", label_x=0.05)
p <- plot_grid(ptmp, get_legend(p3 + theme(legend.position="bottom")), nrow=2, rel_heights = c(2, 0.1))
save_plot(paste0(folder, "IValue.pdf", sep=""), p, base_aspect_ratio = 1.6, base_height=9)

