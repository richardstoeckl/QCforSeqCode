# Copyright Richard Stöckl 2024.
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE or copy at 
# https://www.boost.org/LICENSE_1_0.txt)

library(tidyverse)
library(fs)
library(tinytable)
library(markdown)

# R rounds 0.5 to the nearest even number (e.g. 11.5 -> 12; but 12.5 -> 12), which is not intuitive for most people. 
# This function rounds 0.5 to the nearest higher number. This could influence downstream statistics on the rounded numbers,
# however here, we are only using it for display purposes so I think it is more intuitive.
# the function was implemented from https://stackoverflow.com/questions/12688717/round-up-from-5/12688836#12688836
round2 = function(x, digits) {
    posneg = sign(x)
    z = abs(x)*10^digits
    z = z + 0.5 + sqrt(.Machine$double.eps)
    z = trunc(z)
    z = z/10^digits
    z*posneg
}

# argparse
args<-commandArgs(trailingOnly = TRUE)
checkm2_path <- args[1]
seqkit_stats_path <- args[2]
samtools_coverage_path <- args[3]
gtdbtk_path <- args[4]
decipher_path <- args[5]
tRNAscan_path <- args[6]
outfileReport <- args[7]

# checkm2
checkm2 <- read_tsv(checkm2_path) %>%
            mutate(assembly = Name) %>%
            select(assembly, Completeness,Contamination,Contig_N50)


# seqkit_stats
seqkit_stats_files <- dir_ls(seqkit_stats_path, regexp = "*assembly.tsv", recurse = TRUE)
seqkit_stats <- map_df(seqkit_stats_files,read_tsv, .id = "file") %>% 
                mutate(assembly = str_remove(fs::path_file(fs::as_fs_path(file)),"_assembly.tsv")) %>% 
                select(assembly, num_seqs, max_len)

# samtools coverage
samtools_coverage_files <- dir_ls(samtools_coverage_path, regexp = "*coverage.txt", recurse = TRUE)
samtools_coverage <- map_df(samtools_coverage_files,read_tsv, .id = "file") %>% mutate(assembly = str_remove(fs::path_file(fs::as_fs_path(file)),"_coverage.txt")) %>% select(-file)
samtools_coverage_summary <- samtools_coverage %>% 
                             group_by(assembly) %>% 
                             summarise(mean_depth = mean(meandepth), min_depth = min(meandepth))

# gtdbtk
gtdbtk <- read_tsv(gtdbtk_path) %>% 
    select("Sample ID"=user_genome, "GTDBTK taxonomy"=classification) %>%
    mutate(`GTDBTK taxonomy`=str_remove_all(`GTDBTK taxonomy`, "d__|p__|c__|o__|f__|g__|s__")) %>%
    mutate(`GTDBTK taxonomy`=str_remove_all(`GTDBTK taxonomy`,"_[A-Z]")) %>%
    tidyr::separate("GTDBTK taxonomy", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE)

# decipher
decipher <- read_tsv(decipher_path)

# tRNAscan
tRNAscan_files <- dir_ls(tRNAscan_path, regexp = "*.txt", recurse = TRUE)
tRNAscan <- map_df(tRNAscan_files,read_tsv, .id = "file", skip=3, col_names=c("Sequence Name",
                                                                              "tRNA #",
                                                                              "tRNA begin",
                                                                              "tRNA end",
                                                                              "tRNA type",
                                                                              "Anti codon",
                                                                              "Intron begin",
                                                                              "Intron end",
                                                                              "Inf Score",
                                                                              "Note")) %>%
    mutate(assembly = str_remove(fs::path_file(fs::as_fs_path(file)),".txt")) %>%
    distinct(assembly, `tRNA type`) %>%
    group_by(assembly) %>%
    summarise("tRNA count"=n())


combinedDF <- checkm2 %>% 
    left_join(seqkit_stats, by = join_by(assembly == assembly)) %>% 
    left_join(samtools_coverage_summary, by = join_by(assembly == assembly)) %>%
    left_join(gtdbtk %>% select(`Sample ID`, "GTDB-Tk family"=Family, "GTDB-Tk genus"=Genus), by = join_by(assembly == `Sample ID`)) %>%
    left_join(decipher %>% select(sample, "SILVA family"=family,"SILVA genus"=genus), by = join_by(assembly == sample)) %>%
    left_join(tRNAscan, by = join_by(assembly == assembly)) %>%
    select("Sample ID"=assembly, 
           Completeness, 
           Contamination,
           "Mean coverage depth"=mean_depth, 
           "Minimum coverage depth"=min_depth,
           "GTDB-Tk family",
           "GTDB-Tk genus",
           "SILVA family",
           "SILVA genus",
           "tRNA count",
           "Number of contigs"=num_seqs,
           "N50"=Contig_N50,
           "Length of longest contig"=max_len)


style_table <- function(combinedDF) {
    # convert DF to tinytable
    tab <- tt(combinedDF, tinytable_tt_digits = 1, tinytable_format_num_mark_dec = ".", theme = "bootstrap",
              caption=paste0('Fields marked in <span style="background-color: darkred; color: white;">RED</span> fail to meet the <b>required</b> criteria, ',
                             'Fields marked in <span style="background-color: orange; color: white;">ORANGE</span> fail to meet the <b>recommended</b> criteria ',
                             'as outlined in <a href="https://registry.seqco.de/page/seqcode#data-quality-necessary-for-completion-of-seqcode-registryb" target="_blank">APPENDIX I</a>. ',
                             'Fields marked in <span style="background-color: yellow; color: black;">YELLOW</span> show discrepancies between ',
                             'the whole-genome taxonomic assessment by <a href="https://github.com/Ecogenomics/GTDBTk/" target="_blank">GTDB-Tk</a> and the 16S rRNA gene based taxonomic assessment via <a href="http://www2.decipher.codes/index.html" target="_blank">DECIPHER</a> and the <a href="https://www.arb-silva.de/" target="_blank">SILVA database</a>. ',
                             '<b>Taxonomic assignments are difficult</b>, and I am just a basic automated pipeline, so you should probably <b>look at these samples in detail</b> yourself :).')) %>%
        # add highlight on hover 
        style_tt(bootstrap_class = "table table-hover") %>%
        # and apply some formating for things like decimal numbers, rounding, units, etc.
        format_tt(replace = "-") %>%
        format_tt(j = "Completeness", fn = scales::label_percent(scale = 1)) %>%
        format_tt(j = "Contamination", fn = scales::label_percent(scale=1)) %>%
        format_tt(j = "Mean coverage depth", fn = function(x) paste0(round2(x,digits = 1), "x")) %>%
        format_tt(j = "Minimum coverage depth", fn = function(x) paste0(round2(x,digits = 1), "x")) %>%
        format_tt(j = "N50", fn = function(x) ifelse(x >= 1000000,paste0(round2(x/1000000, digits=2), " Mb"),ifelse(x >= 1000,paste0(round2(x/1000, digits=2), " kb"), x)), digits = 1) %>%
        format_tt(j = "Length of longest contig", fn = function(x) ifelse(x >= 1000000,paste0(round2(x/1000000, digits=2), " Mb"),ifelse(x >= 1000,paste0(round2(x/1000, digits=2), " kb"), x)), digits = 1)
    
    # Apply styles (=background colours) based on conditions
    if (any(combinedDF$Completeness < 90)) {
        tab <- style_tt(tab, i = which(combinedDF$Completeness < 90), j = "Completeness", background = "darkred", color = "white")
    }
    if (any(combinedDF$Contamination > 5)) {
        tab <- style_tt(tab, i = which(combinedDF$Contamination > 5), j = "Contamination", background = "darkred", color = "white")
    }
    if (any(combinedDF$N50 < 25000)) {
        tab <- style_tt(tab, i = which(combinedDF$N50 < 25000), j = "N50", background = "orange", color = "white")
    }
    if (any(combinedDF$`Number of contigs` > 100)) {
        tab <- style_tt(tab, i = which(combinedDF$`Number of contigs` > 100), j = "Number of contigs", background = "orange", color = "white")
    }
    if (any(combinedDF$`Length of longest contig` < 100000)) {
        tab <- style_tt(tab, i = which(combinedDF$`Length of longest contig` < 100000), j = "Length of longest contig", background = "orange", color = "white")
    }
    if (any(combinedDF$`Mean coverage depth` < 10)) {
        tab <- style_tt(tab, i = which(combinedDF$`Mean coverage depth` < 10), j = "Mean coverage depth", background = "darkred", color = "white")
    }
    if (any(combinedDF$`Minimum coverage depth` < 5)) {
        tab <- style_tt(tab, i = which(combinedDF$`Minimum coverage depth` < 5), j = "Minimum coverage depth", background = "darkred", color = "white")
    }
    if (any(combinedDF$`GTDB-Tk family` != combinedDF$`SILVA family`)) {
        tab <- style_tt(tab, i = which(combinedDF$`GTDB-Tk family` != combinedDF$`SILVA family`), j = c("GTDB-Tk family","SILVA family"), background = "yellow", color = "black")
    }
    if (any(combinedDF$`GTDB-Tk genus` != combinedDF$`SILVA genus`)) {
        tab <- style_tt(tab, i = which(combinedDF$`GTDB-Tk genus` != combinedDF$`SILVA genus`), j = c("GTDB-Tk genus","SILVA genus"), background = "yellow", color = "black")
    }
    if (any(combinedDF$`tRNA count` < 18)) {
        tab <- style_tt(tab, i = which(combinedDF$`tRNA count` < 18), j = "tRNA count", background = "orange", color = "white")
    }
    
    # Check for red conditions first
    red_conditions <- combinedDF$Completeness < 90 | combinedDF$Contamination > 5 | combinedDF$`Mean coverage depth` < 10 | combinedDF$`Minimum coverage depth` < 5
    if (any(red_conditions)) {
        tab <- style_tt(tab, i = which(red_conditions), j = 1, background = "darkred", color = "white")
    } else {
        # Check for orange conditions if no red conditions are met
        orange_conditions <- combinedDF$N50 < 25000 | combinedDF$`Number of contigs` > 100 | combinedDF$`Length of longest contig` < 100000
        if (any(orange_conditions)) {
            tab <- style_tt(tab, i = which(orange_conditions), j = 1, background = "orange", color = "white")
        } else {
            # Check for yellow conditions if no red or orange conditions are met
            yellow_conditions <- combinedDF$`GTDB-Tk family` != combinedDF$`SILVA family` | combinedDF$`GTDB-Tk genus` != combinedDF$`SILVA genus`
            if (any(yellow_conditions)) {
                tab <- style_tt(tab, i = which(yellow_conditions), j = 1, background = "yellow", color = "black")
            }
        }
    }
    
    
    return(tab)
}

tab <- style_table(combinedDF)
tab

save_tt(tab, outfileReport, overwrite = TRUE)


