# filter mitochondria and chloroplasts!!!

# package versions
# warnings about phyloseq biom import
# omitted snippets:
# could show snippets of otu counts, sample data, taxon data...



# my.exp <- "15_Hil_CS"
# my.domain <- "Bacteria"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/15_Hil_CS/Bacteria/16S"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "LH/30"



my.exp <- "15_Hil_CS"
my.domain <- "Fungi"
qiime.path <-
  "/home/mark/gitstage/uderica/nordqiime/15_Hil_CS/Fungi/ITS"
my.fact <- "Treatment"
my.ref <- "C/0"
my.other <- "LH/30"



# my.exp <- "15_Hil_HMC"
# my.domain <- "Bacteria"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/15_Hil_HMC/Bacteria/16S"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "LH/30"



# my.exp <- "15_Hil_HMC"
# my.domain <- "Fungi"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/15_Hil_HMC/Fungi/ITS"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "LH/30"



# my.exp <- "Aer_Comp"
# my.domain <- "Fungi"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/Aer_Comp/Aer_Comp_Fungi/raw_data"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "S2-3H-LATE/56"



# my.exp <- "LalStress"
# my.domain <- "Bacteria"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/Lal_Stress/Silva_5101Raw05052017/Bacteria/selected_raw"
# my.fact <- "treatment"
# my.ref <- "C"
# my.other <- "LB500-EARLY"



# my.exp <- "LalStress"
# my.domain <- "Fungi"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/Lal_Stress/Silva_5101Raw05052017/Fungi/raw_data"
# my.fact <- "treatment"
# my.ref <- "C"
# my.other <- "LB500-EARLY"



# my.exp <- "V-HMC"
# my.domain <- "Bacteria"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/V-HMC/Bacteria/raw_data_bacteria"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "LB-AS/90+50h"



# my.exp <- "V-HMC"
# my.domain <- "Fungi"
# qiime.path <-
#   "/home/mark/gitstage/uderica/nordqiime/V-HMC/Fungi/fungi_raw_data"
# my.fact <- "Treatment"
# my.ref <- "C/0"
# my.other <- "LB-AS/90+50h"



diffy.alpha <- 0.05
ord.connect.grps <- TRUE
ord.high.rank <- "Phylum"
ord.low.rank <- "Class"
N.val  <- 10
factor.val <- 100

alpha.div.measures <-
  c("Observed", "Chao1", "Shanon", "Simpson", "InvSimpson")
alpha.div.plot.alab.fontsize <- 6
alpha.point.size <- 5
alpha.point.alpha <- 0.7


prune.chloro <- TRUE
prune.mito <- TRUE


setwd('/home/mark/gitstage/uderica/nordqiime/')



my.timestamp <- format(Sys.time(), format = "%Y%m%d%H%M")

sinkfile <- paste0(my.exp,
                   "_",
                   my.domain,
                   "_",
                   my.timestamp,
                   ".txt")

con <- file(sinkfile)
sink(con, append = TRUE)
sink(con, append = TRUE, type = "message")

pdf(paste0(my.exp,
           "_",
           my.domain,
           "_",
           my.timestamp,
           ".pdf"))

source("uderica_fxns.R")


# # first load data
# phylo <- uderica.load(the.exp = my.exp, the.domain = my.domain)

phylo <-
  nordqiime.load(qiime.path,
                 my.domain)

x <- as.data.frame(phylo@tax_table@.Data)

table(x$Class == "Chloroplast")
table(x$Family == "mitochondria")
chloro.flag <- x$Class == "Chloroplast"
mito.flag <- x$Family == "mitochondria"

if (prune.chloro) {
  phylo <- prune_taxa(!chloro.flag, phylo)
}

if (prune.mito) {
  phylo <- prune_taxa(!mito.flag, phylo)
}



###   ###   ###
# make me a function
otu.dat <- phylo@otu_table
dim(otu.dat)

# could show snippets of otu counts
# otu.tab.excerpt <- otu.dat[order(rownames(otu.dat)),]
# otu.tab.excerpt <- otu.tab.excerpt[,order(as.numeric(colnames(otu.dat)))]
# pander(otu.tab.excerpt[1:9,1:9])

unlisted.otu.counts <- unlist(as.matrix.data.frame(otu.dat))

old.mar <- par("mar")
par("mar" = c(5, 4, 3, 2))

hist(
  log10(unlisted.otu.counts),
  breaks = 100,
  main = paste0(my.exp, " ", my.domain, "\nHistogram of reads per sample/OTU"),
  xlab = paste0(my.exp, " ", my.domain, " log10(read counts)")
)
par("mar" = old.mar)
###   ###   ###


###   ###   ###
# make me a function
map.dat <- as.data.frame(phylo@sam_data@.Data)
names(map.dat) <- phylo@sam_data@names
mdncol <- ncol(map.dat)
# assumes a particular map.txt structure (has been consistent so far)
exp.design.matrix <- map.dat[, 5:(mdncol - 1)]

print(pander(table(exp.design.matrix)))
###   ###   ###

###   ###   ###
# make me a function
# make a copy of teh taxonomy data
tax.dat <- as.data.frame(phylo@tax_table@.Data)

# make it a matrix... so we can easily test individual values
temp <- as.matrix(tax.dat)

# make the cells NA if the value is unidentified OR starts with uncultured
temp[grepl(pattern = "^uncultured", temp)] <- NA
temp[temp == "unidentified"] <- NA

# put TRUE (1) in each cell if it's not NA
temp <- !is.na(temp)

# sum the TRUEs in each row
# this actually counts the number of determined ranks
# K F C O F G S
# 1 2 3 4 5 6 7
# I figured that, if KFCOF were known, then GS wouldn't be, and the sum would be 5
# but with fungi, in internal rank could be undetermined, followed by a known rank
# KFCO G

temp2 <- tax.dat[, !grepl(pattern = "^Rank", colnames(tax.dat))]
temp2$taxonomic.specificity

taxonomic.specificity <- rowSums(temp)
tax.dat$taxonomic.specificity <- taxonomic.specificity

temp2$taxonomic.specificity <- taxonomic.specificity

max.specificity <- max(tax.dat$taxonomic.specificity)
temp.frame <-
  cbind.data.frame((1:max.specificity), names(tax.dat)[1:max.specificity])
names(temp.frame) <- c("rank", "label")

max.count <- hist(tax.dat$taxonomic.specificity, plot = FALSE)
max.count <- max.count$counts
max.count <-  max(max.count)

par("mar" = c(5, 4, 3, 2))

hist(
  tax.dat$taxonomic.specificity,
  labels = as.character(temp.frame$label),
  breaks = 0:max.specificity,
  main = paste0(
    my.exp,
    " ",
    my.domain,
    "\nHistogram of INFERRED most specific rank"
  ),
  xlab = 'Number of non "uncultured" and "unidentified" phylogenetic levels',
  ylim = c(0, (max.count + 20))
)

par("mar" =  old.mar)
###   ###   ###



# optionally filter
phylo.top.N <- top.N.filter(the.phylo = phylo, the.N = N.val)
phylo.k.of.A <-
  k.of.A.filter(the.phylo = phylo, k.factor = factor.val)

# barplot
uderica.barplot(
  the.phylo = phylo.top.N,
  the.factor = my.fact,
  the.rank = "Family",
  the.filter.meth = "top n"
)
uderica.barplot(
  the.phylo = phylo.top.N,
  the.factor = my.fact,
  the.rank = "Genus",
  the.filter.meth = "k of a"
)
uderica.barplot(
  the.phylo = phylo.k.of.A,
  the.factor = my.fact,
  the.rank = "Family",
  the.filter.meth = "top n"
)
uderica.barplot(
  the.phylo = phylo.k.of.A,
  the.factor = my.fact,
  the.rank = "Genus",
  the.filter.meth = "k of a"
)


# alpha/richness
uderica.alpha(the.phylo = phylo,
              the.factor = my.fact,
              the.measures = alpha.div.measures)

# ordination/"CoA"
ord_meths  <- c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
uderica.ordinate(
  the.phylo = phylo,
  the.meths = ord_meths,
  the.dist = "bray",
  the.fact = my.fact,
  connect = ord.connect.grps
)

# other distances, some require addirtional arguments
# bray gower (daisy?) and jsd seem to work
# "unifrac"
# Original (unweighted) UniFrac distance, UniFrac
# "wunifrac"
# weighted-UniFrac distance, UniFrac
# "dpcoa"
# sample-wise distance used in Double Principle Coordinate Analysis, DPCoA
# "jsd"


# heatmap, most distant from control
constrast.helper(the.phylo = phylo,
                 the.fact = my.fact,
                 the.ref = my.ref)

# differentail analysis
uderica.diffy(
  the.phylo = phylo,
  the.factor = my.fact,
  the.contrast.num = my.ref,
  the.contrast.denom =  my.other,
  the.alpha = diffy.alpha,
  the.high.rank = ord.high.rank,
  the.low.rank = ord.low.rank,
  the.y.axis = "log2FoldChange"
)


dev.off()

sink()
sink(type = "message")
