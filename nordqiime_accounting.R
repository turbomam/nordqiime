library(plyr)
library(readxl)
library(readr)

# SO MUCH REDUNDANT CODE!
# assumptions:  only one study or domain chunk will appear in each path
# study names may vary a little here and there,
# but will always have the same letters and numbers in the same order
# no insertions, deletions or rearrangements of letters and numbers
# only different punctuation

cat.and.label.biom.sum <-
  function(root.path) {
    desired.files <-
      list.files(
        path = root.path,
        pattern = "biom_table_summary.txt",
        full.names = TRUE,
        recursive = TRUE
      )
    
    the.agg <- lapply(desired.files, function(one.desired) {
      # print(one.desired)
      desired.chunks <- strsplit(one.desired, "/")
      desired.rhs <-
        desired.chunks[[1]][length(desired.chunks[[1]])]
      # print(c(desired.rhs, desired.name))
      exact.match  <-  desired.rhs == "biom_table_summary.txt"
      if (exact.match) {
        curent.read <-   read_delim(
          one.desired,
          ":",
          escape_double = FALSE,
          col_names = FALSE,
          trim_ws = TRUE
        )
        
      }
      
      curent.read$source.file <- one.desired
      return(curent.read)
      
    })
    the.agg <- do.call(rbind.fill, the.agg)
    names(the.agg) <- make.names(names(the.agg))
    return(the.agg)
  }


cat.and.label.txt.or.xlsx <-
  function(desired.name,
           root.path,
           mismatch.allowed) {
    desired.files <-
      list.files(
        path = root.path,
        pattern = desired.name,
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = mismatch.allowed
      )
    
    the.agg <- lapply(desired.files, function(one.desired) {
      # print(one.desired)
      desired.chunks <- strsplit(one.desired, "/")
      desired.rhs <-
        desired.chunks[[1]][length(desired.chunks[[1]])]
      # print(c(desired.rhs, desired.name))
      exact.match  <-  desired.rhs == desired.name
      if (exact.match || mismatch.allowed) {
        # print("identical")
        name.chunks <- strsplit(desired.rhs, "\\.")
        # print(name.chunks)
        current.extension <-
          name.chunks[[1]][length(name.chunks[[1]])]
        # print(current.extension)
        if (current.extension == "txt") {
          curent.read <-
            read.delim(one.desired)
        } else if (current.extension == "xlsx") {
          curent.read <- read_excel(one.desired)
          
        }
        
        curent.read$source.file <- one.desired
        return(curent.read)
      } else {
        # print("not quite")
      }
      
    })
    the.agg <- do.call(rbind.fill, the.agg)
    names(the.agg) <- make.names(names(the.agg))
    return(the.agg)
  }

multiqc.agg <-
  cat.and.label.txt.or.xlsx("multiqc_fastqc.txt",
                            "/home/mark/gitstage/uderica/nordqiime/",
                            FALSE)

multiqc.agg$Sample <- as.character(multiqc.agg$Sample)

study.opts <-
  list.dirs(
    "/home/mark/gitstage/uderica/nordqiime",
    recursive = FALSE,
    full.names = FALSE
  )

study.opts <- setdiff(study.opts, "FastQC")

lapply(study.opts, function(one.study) {
  study.flag <-
    grepl(pattern = one.study,
          x = multiqc.agg$source.file,
          ignore.case = TRUE)
  multiqc.agg$study[study.flag] <<- one.study
})

multiqc.agg$study.letters <- tolower(multiqc.agg$study)
multiqc.agg$study.letters <- make.names(multiqc.agg$study.letters)
multiqc.agg$study.letters <-
  gsub(pattern = "\\.+", "", multiqc.agg$study.letters)
multiqc.agg$study.letters <-
  gsub(pattern = "\\_+", "", multiqc.agg$study.letters)
multiqc.agg$study.letters <-
  sub(pattern = "^X", "", multiqc.agg$study.letters)

domain.opts <- c("bacteria", "fungi")

lapply(domain.opts, function(one.domain) {
  domain.flag <-
    grepl(pattern = one.domain,
          x = multiqc.agg$source.file,
          ignore.case = TRUE)
  multiqc.agg$domain[domain.flag] <<- one.domain
})

multiqc.agg$ID <- sub(
  pattern = "\\.extendedFrags$",
  replacement = "",
  x = multiqc.agg$Sample,
  ignore.case = FALSE
)

multiqc.agg$ID[!grepl(pattern = "\\.extendedFrags$", multiqc.agg$Sample)] <-
  NA



flashed.flag <- grepl(pattern = "extended", x = multiqc.agg$Sample)

initial.multiqc.agg <- multiqc.agg[!flashed.flag, ]
flashed.multiqc.agg <- multiqc.agg[flashed.flag, ]

# map.agg <-
#   cat.and.label.txt.or.xlsx("map.txt",
#                             "/home/mark/gitstage/uderica/nordqiime/", FALSE)

xlsx.agg <-
  cat.and.label.txt.or.xlsx("file path.xlsx",
                            "/home/mark/gitstage/uderica/nordqiime/",
                            TRUE)

xlsx.agg$study.letters <- tolower(xlsx.agg$study)
xlsx.agg$study.letters <- make.names(xlsx.agg$study.letters)
xlsx.agg$study.letters <-
  gsub(pattern = "\\.+", "", xlsx.agg$study.letters)
xlsx.agg$study.letters <-
  sub(pattern = "^X", "", xlsx.agg$study.letters)

lapply(domain.opts, function(one.domain) {
  domain.flag <-
    grepl(pattern = one.domain,
          x = xlsx.agg$source.file,
          ignore.case = TRUE)
  xlsx.agg$domain[domain.flag] <<- one.domain
})

xlsx.agg$key <-
  paste(xlsx.agg$study.letters, xlsx.agg$domain, xlsx.agg$ID, sep = ".")

xlsx.agg$File.R1.base <-
  sub(pattern = "\\..*$", "", x =  xlsx.agg$File.R1)
xlsx.agg$File.R2.base <-
  sub(pattern = "\\..*$", "", x =  xlsx.agg$File.R2)

# unique(xlsx.agg$source.file)

# temp <- table(xlsx.agg$ID)
# temp[temp != 2]
# names(temp[temp != 2])
# sort(names(temp[temp != 2]))
#
# sort(unique(xlsx.agg$source.file))



merged.1 <- merge(
  x = initial.multiqc.agg,
  y = xlsx.agg,
  by.x = "Sample",
  by.y = "File.R1.base",
  all = FALSE
)

merged.2 <- merge(
  x = initial.multiqc.agg,
  y = xlsx.agg,
  by.x = "Sample",
  by.y = "File.R2.base",
  all = FALSE
)

merged.combo <- rbind.fill(merged.1, merged.2)

merged.combo$source.file[!is.na(merged.combo$source.file.x)] <-
  merged.combo$source.file.x[!is.na(merged.combo$source.file.x)]

merged.combo$source.file[!is.na(merged.combo$source.file.y)] <-
  merged.combo$source.file.y[!is.na(merged.combo$source.file.y)]

# collapse file name bases?  still split between R1 and R2

initial.numeric <-
  merged.combo[, c(#    "Sequences.flagged.as.poor.quality",
    "avg_sequence_length",
    "Total.Sequences",
    "key")]

initial.numeric <-
  aggregate(
    initial.numeric,
    by = list(initial.numeric$key),
    FUN = mean,
    na.rm = TRUE
  )

initial.numeric$key <- initial.numeric$Group.1
keepers <- setdiff(names(initial.numeric), "Group.1")
initial.numeric <- initial.numeric[, keepers]

flashed.multiqc.agg$key <-
  paste(
    flashed.multiqc.agg$study.letters,
    flashed.multiqc.agg$domain,
    flashed.multiqc.agg$ID,
    sep = "."
  )

flashed.numeric <-
  flashed.multiqc.agg[, c(#    "Sequences.flagged.as.poor.quality",
    "avg_sequence_length",
    "Total.Sequences",
    "key")]

flash.result <- merge(
  x = initial.numeric,
  y = flashed.numeric,
  by = "key",
  suffixes = c(".initial", ".flashed")
)

flash.result$flash.pct.calc <- 1 - (
  (
    flash.result$Total.Sequences.initial - flash.result$Total.Sequences.flashed
  ) / flash.result$Total.Sequences.initial
)

biom.summaries <-
  cat.and.label.biom.sum("/home/mark/gitstage/uderica/nordqiime/")

agg.flag <- biom.summaries$X1 %in% c(
  "Counts/sample detail",
  "Counts/sample summary",
  "Max",
  "Mean",
  "Median",
  "Min",
  "Num observations",
  "Num samples",
  "Observation Metadata Categories",
  "Sample Metadata Categories",
  "Std. dev.",
  "Table density (fraction of non-zero values)",
  "Total count"
)

sample.biom.summaries <- biom.summaries[!agg.flag, ]

lapply(study.opts, function(one.study) {
  study.flag <-
    grepl(pattern = one.study,
          x = sample.biom.summaries$source.file,
          ignore.case = TRUE)
  sample.biom.summaries$study[study.flag] <<- one.study
})

sample.biom.summaries$study.letters <-
  tolower(sample.biom.summaries$study)
sample.biom.summaries$study.letters <-
  make.names(sample.biom.summaries$study.letters)
sample.biom.summaries$study.letters <-
  gsub(pattern = "\\.+", "", sample.biom.summaries$study.letters)
sample.biom.summaries$study.letters <-
  gsub(pattern = "\\_+", "", sample.biom.summaries$study.letters)
sample.biom.summaries$study.letters <-
  sub(pattern = "^X", "", sample.biom.summaries$study.letters)

lapply(domain.opts, function(one.domain) {
  domain.flag <-
    grepl(pattern = one.domain,
          x = sample.biom.summaries$source.file,
          ignore.case = TRUE)
  sample.biom.summaries$domain[domain.flag] <<- one.domain
})

sample.biom.summaries$key <-
  paste(
    sample.biom.summaries$study.letters,
    sample.biom.summaries$domain,
    sample.biom.summaries$X1,
    sep = "."
  )

# really should get the summary by running XYZ command on LNM bio files
# (different between fungi adn bacteria)
# these cordiv files don't report samples that were rejected due to low counts
sample.biom.summaries <- sample.biom.summaries[, c("key", "X2")]
names(sample.biom.summaries) <- c("key", "Total.Sequences.picked")

flash.result <- merge(x = flash.result,
                      y = sample.biom.summaries,
                      by = "key",
                      all = TRUE)

flash.result$Total.Sequences.picked <-
  as.numeric(flash.result$Total.Sequences.picked)

flash.result$pick.pct.calc <- 1 - (
  (
    flash.result$Total.Sequences.flashed - flash.result$Total.Sequences.picked
  ) / flash.result$Total.Sequences.flashed
)

temp <- strsplit(flash.result$key, "\\.")
temp <- do.call(rbind.data.frame, temp)
names(temp) <- c("study", "domain", "ID")

flash.result <- cbind.data.frame(temp, flash.result)

# spot checked the (hard-to-parse) flash output... their reported "extended" percentages
# match my calculated percentages

# now for blasting the notcombined sequences...
# sed -n '1~4s/^@/>/p;2~4p' 191.notCombined_1.fastq > 191.notCombined_1.fasta

# vhmc.bacteria.1140
# sed -n '1~4s/^@/>/p;2~4p' /home/mark/gitstage/uderica/nordqiime/V-HMC/Bacteria/raw_data_bacteria/extendedFrags/1140.notCombined_1.fastq > /home/mark/gitstage/uderica/nordqiime/V-HMC/Bacteria/raw_data_bacteria/extendedFrags/1140.notCombined_1.fasta
# sed -n '1~4s/^@/>/p;2~4p' /home/mark/gitstage/uderica/nordqiime/V-HMC/Bacteria/raw_data_bacteria/extendedFrags/1140.notCombined_2.fastq > /home/mark/gitstage/uderica/nordqiime/V-HMC/Bacteria/raw_data_bacteria/extendedFrags/1140.notCombined_2.fasta
# cat 1140.notCombined_1.fasta 1140.notCombined_2.fasta > 1140.notCombined_both.fasta
