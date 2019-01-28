library(polyester)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)

fasta = readDNAStringSet(args[1])
# renaming transcript names
oldnames <- strsplit(names(fasta), " ")
newnames <- lapply(oldnames, function(x) x[1])
newnames <- unlist(newnames)
names(fasta) <- newnames

readlen = 100
readspertx = round(1 * width(fasta) / readlen)
fold_changes = matrix(rep(1, 2*length(fasta)), nrow=length(fasta))

theoretical <- read.table(args[2], header=FALSE, sep="\t")
colnames(theoretical) <- c("Name", "Length", "CopyNumber", "NumReads")

theoretical <- theoretical[match(names(fasta), theoretical$Name), ]
for (i in 1:length(fasta)) {
	readspertx[i] = max(round(theoretical[i, "NumReads"]), 1)
}

#simulate_experiment(args[1], reads_per_transcript=readspertx, num_reps=c(1,1), fold_changes=fold_changes, error_model='illumina5', 
#	bias='cdnaf', outdir=args[3])
simulate_experiment(args[1], reads_per_transcript=readspertx, num_reps=c(1,1), fold_changes=fold_changes, error_model='illumina5', outdir=args[3])

