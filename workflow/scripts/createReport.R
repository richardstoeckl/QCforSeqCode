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
outfileReportHTML <- args[7]
outfileReportTSV <- args[8]
samples <- read_csv(args[9], col_types = cols(sampleID = col_character()))$sampleID

# checkm2
checkm2 <- read_tsv(checkm2_path) %>%
            mutate(assembly = Name) %>%
            select(assembly, Completeness,Contamination,Contig_N50)


# seqkit_stats
seqkit_stats_files <- file.path(seqkit_stats_path, paste0(samples, "_assembly.tsv")) %>% set_names()
seqkit_stats <- map_df(seqkit_stats_files,read_tsv, .id = "file") %>% 
                mutate(assembly = str_remove(fs::path_file(fs::as_fs_path(file)),"_assembly.tsv")) %>% 
                select(assembly, num_seqs, max_len)

# samtools coverage
samtools_coverage_files <- file.path(samtools_coverage_path, paste0(samples, "_coverage.txt")) %>% set_names()
samtools_coverage <- map_df(samtools_coverage_files,read_tsv, .id = "file") %>% mutate(assembly = str_remove(fs::path_file(fs::as_fs_path(file)),"_coverage.txt")) %>% select(-file)
samtools_coverage_summary <- samtools_coverage %>% 
                             mutate(contig_length = endpos - startpos + 1) %>%
                             group_by(assembly) %>% 
                             summarise(mean_depth = weighted.mean(meandepth, contig_length),
                                       min_contig_mean_depth = min(meandepth), .groups = "drop")

# gtdbtk
gtdbtk <- read_tsv(gtdbtk_path) %>% 
    select("Sample ID"=user_genome, "GTDBTK taxonomy"=classification) %>%
    mutate(`GTDBTK taxonomy`=str_remove_all(`GTDBTK taxonomy`, "d__|p__|c__|o__|f__|g__|s__")) %>%
    mutate(`GTDBTK taxonomy`=str_remove_all(`GTDBTK taxonomy`,"_[A-Z]")) %>%
    tidyr::separate("GTDBTK taxonomy", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE)

# decipher
decipher <- read_tsv(decipher_path)

# tRNAscan
tRNAscan_files <- file.path(tRNAscan_path, paste0(samples, ".txt")) %>% set_names()
tRNAscan <- map_df(tRNAscan_files,read_tsv, .id = "file", skip=3, col_types=cols(.default = col_character()), col_names=c("Sequence Name",
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


combinedDF <- tibble(assembly = samples) %>%
    left_join(checkm2, by = join_by(assembly == assembly)) %>%
    left_join(seqkit_stats, by = join_by(assembly == assembly)) %>% 
    left_join(samtools_coverage_summary, by = join_by(assembly == assembly)) %>%
    left_join(gtdbtk %>% select(`Sample ID`, "GTDB-Tk family"=Family, "GTDB-Tk genus"=Genus), by = join_by(assembly == `Sample ID`)) %>%
    left_join(decipher %>% select(sample, "SILVA family"=family,"SILVA genus"=genus), by = join_by(assembly == sample)) %>%
    left_join(tRNAscan, by = join_by(assembly == assembly)) %>%
    mutate(`tRNA count` = replace_na(`tRNA count`, 0L)) %>%
    mutate("GTDB-Tk family"=ifelse(is.na(`GTDB-Tk family`), "UNCLASSIFIED", `GTDB-Tk family`),
        "GTDB-Tk genus"=ifelse(is.na(`GTDB-Tk genus`), "UNCLASSIFIED", `GTDB-Tk genus`),
        "SILVA family"=ifelse(is.na(`SILVA family`), "UNCLASSIFIED", `SILVA family`),
        "SILVA genus"=ifelse(is.na(`SILVA genus`), "UNCLASSIFIED", `SILVA genus`)) %>%
    select("Sample ID"=assembly, 
           Completeness, 
           Contamination,
           "Mean coverage depth"=mean_depth, 
           "Lowest contig mean depth"=min_contig_mean_depth,
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
    tab <- tt(combinedDF,
              caption=paste0('Fields marked in <span style="background-color: darkred; color: white;">RED</span> fail to meet the <b>required</b> criteria, ',
                             'Fields marked in <span style="background-color: orange; color: white;">ORANGE</span> fail to meet the <b>recommended</b> criteria ',
                             'as outlined in <a href="https://registry.seqco.de/page/seqcode#data-quality-necessary-for-completion-of-seqcode-registryb" target="_blank">APPENDIX I</a>. ',
                             'Fields marked in <span style="background-color: yellow; color: black;">YELLOW</span> show discrepancies between ',
                             'the whole-genome taxonomic assessment by <a href="https://github.com/Ecogenomics/GTDBTk/" target="_blank">GTDB-Tk</a> and the 16S rRNA gene based taxonomic assessment via <a href="http://www2.decipher.codes/index.html" target="_blank">DECIPHER</a> and the <a href="https://www.arb-silva.de/" target="_blank">SILVA database</a>. ',
                             '<b>Taxonomic assignments are difficult</b>, and I am just a basic automated pipeline, so you should probably <b>look at these samples in detail</b> yourself :).')) %>%
        theme_html(
            engine = "bootstrap",
            class = "table table-hover",
            portable = TRUE
        ) %>%
        # and apply some formating for things like decimal numbers, rounding, units, etc.
        format_tt(replace = "-") %>%
        format_tt(j = "Completeness", fn = scales::label_percent(scale = 1)) %>%
        format_tt(j = "Contamination", fn = scales::label_percent(scale=1)) %>%
        format_tt(j = "Mean coverage depth", fn = function(x) paste0(round2(x,digits = 1), "x")) %>%
        format_tt(j = "Lowest contig mean depth", fn = function(x) paste0(round2(x,digits = 1), "x")) %>%
        format_tt(j = "N50", fn = function(x) ifelse(x >= 1000000,paste0(round2(x/1000000, digits=2), " Mb"),ifelse(x >= 1000,paste0(round2(x/1000, digits=2), " kb"), x)), digits = 1) %>%
        format_tt(j = "Length of longest contig", fn = function(x) ifelse(x >= 1000000,paste0(round2(x/1000000, digits=2), " Mb"),ifelse(x >= 1000,paste0(round2(x/1000, digits=2), " kb"), x)), digits = 1)
    
    # Apply styles (=background colours) based on conditions
    if (any(combinedDF$Completeness < 90, na.rm = TRUE)) {
        tab <- style_tt(tab, i = which(combinedDF$Completeness < 90), j = "Completeness", background = "darkred", color = "white")
    }
    if (any(combinedDF$Contamination > 5, na.rm = TRUE)) {
        tab <- style_tt(tab, i = which(combinedDF$Contamination > 5), j = "Contamination", background = "darkred", color = "white")
    }
    if (any(combinedDF$N50 < 25000, na.rm = TRUE)) {
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
    if (any(combinedDF$`Lowest contig mean depth` < 5)) {
        tab <- style_tt(tab, i = which(combinedDF$`Lowest contig mean depth` < 5), j = "Lowest contig mean depth", background = "darkred", color = "white")
    }
    if (any(combinedDF$`GTDB-Tk family` != combinedDF$`SILVA family`)) {
        tab <- style_tt(tab, i = which(combinedDF$`GTDB-Tk family` != combinedDF$`SILVA family`), j = c("GTDB-Tk family","SILVA family"), background = "yellow", color = "black")
    }
    if (any(combinedDF$`GTDB-Tk genus` != combinedDF$`SILVA genus`)) {
        tab <- style_tt(tab, i = which(combinedDF$`GTDB-Tk genus` != combinedDF$`SILVA genus`), j = c("GTDB-Tk genus","SILVA genus"), background = "yellow", color = "black")
    }
    if (any(combinedDF$`GTDB-Tk genus` == "UNCLASSIFIED")) {
        tab <- style_tt(tab, i = which(combinedDF$`GTDB-Tk genus` == "UNCLASSIFIED"), j = "GTDB-Tk genus", color = "orange")
    }
    if (any(combinedDF$`GTDB-Tk family` == "UNCLASSIFIED")) {
        tab <- style_tt(tab, i = which(combinedDF$`GTDB-Tk family` == "UNCLASSIFIED"), j = "GTDB-Tk family", color = "orange")
    }
    if (any(combinedDF$`SILVA genus` == "UNCLASSIFIED")) {
        tab <- style_tt(tab, i = which(combinedDF$`SILVA genus` == "UNCLASSIFIED"), j = "SILVA genus", color = "orange")
    }
    if (any(combinedDF$`SILVA family` == "UNCLASSIFIED")) {
        tab <- style_tt(tab, i = which(combinedDF$`SILVA family` == "UNCLASSIFIED"), j = "SILVA family", color = "orange")
    }
    if (any(combinedDF$`tRNA count` < 18)) {
        tab <- style_tt(tab, i = which(combinedDF$`tRNA count` < 18), j = "tRNA count", background = "orange", color = "white")
    }
    
    # Check for red conditions first
    red_conditions <- combinedDF$Completeness < 90 | combinedDF$Contamination > 5 | combinedDF$`Mean coverage depth` < 10 | combinedDF$`Lowest contig mean depth` < 5
    if (any(red_conditions, na.rm = TRUE)) {
        tab <- style_tt(tab, i = which(red_conditions), j = 1, background = "darkred", color = "white")
    } else {
        # Check for orange conditions if no red conditions are met
        orange_conditions <- combinedDF$N50 < 25000 | combinedDF$`Number of contigs` > 100 | combinedDF$`Length of longest contig` < 100000
        if (any(orange_conditions, na.rm = TRUE)) {
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

save_tt(tab, outfileReportHTML, overwrite = TRUE)

write_tsv(combinedDF, outfileReportTSV)


