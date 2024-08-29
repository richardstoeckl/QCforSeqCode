# Copyright Richard Stöckl 2024.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at 
# https://www.boost.org/LICENSE_1_0.txt)


library(tidyverse)
library(fs)
library(DECIPHER)


# argparse
args <- commandArgs(trailingOnly = TRUE)
trainingData <- args[1]
fastaPath <- args[2]
threads <- args[3]
outPath <- args[4]

# load trained data
load(trainingData)

# load fasta file
fasta <- readDNAStringSet(fastaPath)

# assign taxonomic classifiction
taxonomy <- IdTaxa(fasta, training = trainingSet, processors = as.numeric(threads), type= "extended")

# Initialize an empty list to store the results
results_list <- list()

# Loop through the taxonomy data
for (i in seq_along(taxonomy)) {
    assignment <- taxonomy[[i]]$taxon
    conf <- taxonomy[[i]]$confidence
    rank <- taxonomy[[i]]$rank
    tibble <- tibble(assignment = assignment, confidence = min(conf), rank = rank)
    tibble_pivot <- tibble %>%
        pivot_wider(names_from = rank, values_from = assignment) %>% 
        mutate(sample = names(taxonomy)[i])
    # Append the result to the list
    results_list[[i]] <- tibble_pivot
}

# Combine all tibbles into a single dataframe
final_dataframe <- bind_rows(results_list)

write_tsv(final_dataframe, outPath)





