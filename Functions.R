library(tidyverse)
library(data.table)
library(stringr)
library(readxl)
library(readr)
library(combinat)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(GGally)
library(forcats)
library(cowplot)
library(bio3d)
library(ggseqlogo)
library(ggsignif)
library(ggh4x)

############################## func 1 ###################################
# to check mutational frequency in files Rep0/1/2 P0/4, where there is a spike ~ 570 nt
mutfreq <- function(filename){
  mut_freq <- read_csv(paste0(filedir,filename,"_codoncounts.csv")) %>% 
    filter(site > 550 & site <600) %>%
    left_join(amplicon, by = "site") %>% 
    pivot_longer(cols = -c(site, wildtype, amplicon),
                 names_to = "codon", values_to = "count") %>% 
    mutate(type = case_when(wildtype == codon ~ "wildtype",
                            T ~ "variant")) %>% 
    group_by(site, type) %>% 
    summarize(count = sum(count)) %>% 
    ungroup() %>% 
    group_by(site) %>% 
    spread(key = type, value = count) %>% 
    mutate(mut.freq = variant/(variant + wildtype))
  
  p <- ggplot(mut_freq, aes(x=site, y=mut.freq)) +
    geom_line() +
    ggtitle(paste("Mutational frequency \n of", filename)) +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  
  return(p)
}

############################## func 2 ###################################
# function to calculate how many unique codons are there in site 577 & 578
# to check if there were more codons being made because of sample prep
mut_codon_count <- function(filename){
  mut_codoncount <- read_csv(paste0(filedir,filename,"_codoncounts.csv")) %>% 
    filter(site >= 577 & site <= 578) %>%
    left_join(amplicon, by = "site") %>% 
    pivot_longer(cols = -c(site, wildtype, amplicon),
                 names_to = "codon", values_to = "count") %>% 
    mutate(type = case_when(wildtype == codon ~ "wildtype",
                            T ~ "variant")) %>% 
    group_by(site) %>% 
    summarize(total.codon = sum(count != 0))
  
  return(mut_codoncount)
}

########################### helper df 1 #################################
codons <- tibble(
  codon = c("TTT","TTC",
            "TTA","TTG","CTT","CTC","CTA","CTG",
            "ATT","ATC","ATA",
            "ATG",
            "GTT","GTC","GTA","GTG",
            "TCT","TCC","TCA","TCG","AGT","AGC",
            "CCT","CCC","CCA","CCG",
            "ACT","ACC","ACA","ACG",
            "GCT","GCC","GCA","GCG",
            "TAT","TAC",
            "TAA","TAG","TGA",
            "CAT","CAC",
            "CAA","CAG",
            "AAT","AAC",
            "AAA","AAG",
            "GAT","GAC",
            "GAA","GAG",
            "TGT","TGC",
            "TGG",
            "CGT","CGC","CGA","CGG","AGA","AGG",
            "GGT","GGC","GGA","GGG"),
  aa = c("Phe","Phe",
         "Leu","Leu","Leu","Leu","Leu","Leu",
         "Ile","Ile","Ile",
         "Met",
         "Val","Val","Val","Val",
         "Ser","Ser","Ser","Ser","Ser","Ser",
         "Pro","Pro","Pro","Pro",
         "Thr","Thr","Thr","Thr",
         "Ala","Ala","Ala","Ala",
         "Tyr","Tyr",
         "***","***","***",
         "His","His",
         "Gln","Gln",
         "Asn","Asn",
         "Lys","Lys",
         "Asp","Asp",
         "Glu","Glu",
         "Cys","Cys",
         "Trp",
         "Arg","Arg","Arg","Arg","Arg","Arg",
         "Gly","Gly","Gly","Gly")
)

############################## func 3 ###################################
assign_amplicon <- function(site_tibble){
  t <- mutate(site_tibble, amplicon = case_when(
    site >= 1 & site <= 93 ~ "1",
    site >= 94 & site <= 189 ~ "2",
    site >= 190 & site <= 290 ~ "3",
    site >= 291 & site <= 400 ~ "4",
    site >= 401 & site <= 503 ~ "5",
    site >= 504 & site <= 592 ~ "6",
    site >= 593 & site <= 686 ~ "7",
    site >= 687 & site <= 758 ~ "8",
    TRUE ~ "invalid"))
  return(t)
}

########################### helper df 2 #################################
amplicon <- assign_amplicon(tibble(site = c(1:758)))

############################## func 4 ###################################
# similar to "read_counts", but modified specifically for calculating diversity decrease
# through passages.
get_counts <- function(filename){
  counts <- read_csv(paste0(filedir,filename, "_codoncounts.csv")) %>%
    left_join(amplicon, by = "site") %>% 
    #translate by joining to codons dataframe
    pivot_longer(cols = -c(site, wildtype, amplicon),
                 names_to = "codon", values_to = "count") %>%
    left_join(codons, by = "codon") %>%
    left_join(codons, by = c("wildtype" = "codon")) %>%
    #assign mutation types
    mutate(mutation = paste0(aa.y, site, aa.x, sep = ""),
           mutation.type = case_when(
             aa.x == aa.y & codon == wildtype ~ "wildtype",
             aa.x == aa.y ~ "silent",
             aa.x == "***" ~ "nonsense",
             TRUE ~ "missense"
           )) %>% 
    filter(site<758)
  
  counts$ending <- substr(gsub('\\d','', counts$codon), 3, 3)
  
  remainingcounts <- counts %>% 
    filter(mutation.type == "wildtype" | ending == "C" | ending == "G") %>%
    group_by(site, amplicon, codon) %>%
    summarize(count = sum(count)) %>% 
    filter(count > 0)
  
  remainingaa<- counts %>% 
    # filter(mutation.type == "wildtype" | ending == "C" | ending == "G") %>% 
    group_by(site, amplicon, mutation) %>%
    summarize(count = sum(count)) %>% 
    filter(count>0)
  
  file_counts <- c(dim(remainingcounts)[1], dim(remainingaa)[1])
  
  return(file_counts)
}

############################## func 5 ###################################
# used for amino acid preference calculation
# examples input: "Rep0Pla"
read_counts <- function(filename){
  counts <- read_csv(paste0(filedir, filename, "_codoncounts.csv")) %>%
    left_join(amplicon, by = "site") %>% 
    #translate by joining to codons dataframe
    pivot_longer(cols = -c(site, wildtype, amplicon),
                 names_to = "codon", values_to = "count") %>%
    left_join(codons, by = "codon") %>%
    left_join(codons, by = c("wildtype" = "codon")) %>%
    #assign mutation types
    mutate(mutation = paste0(aa.y, site, aa.x, sep = ""),
           mutation.type = case_when(
             aa.x == aa.y & codon == wildtype ~ "wildtype",
             aa.x == aa.y ~ "silent",
             aa.x == "***" ~ "nonsense",
             TRUE ~ "missense"
           )) %>%
    #collapse mutations
    group_by(site, amplicon, mutation, mutation.type) %>%
    summarize(count = sum(count))
  
  return(counts)
}

############################## func 6 ###################################
# used in amino acid decrease section rather than amino acid preferences section
# compute enrichment by replicate (instead of pre-, wt-, post-passage)
# compute enrichment pla <- p0, pla <- p1, pla <- p4
compute_enrichment_by_replicate <- function(pla, WT_pla, p0, p1, p4){
  Mut_plasmid <- read_counts(pla)
  WT_plasmid <- read_counts(WT_pla)
  P0 <- read_counts(p0)
  P1 <- read_counts(p1)
  P4 <- read_counts(p4)
  
  #filter for counts >10 in plasmid libary
  Input <- Mut_plasmid %>%
    filter(count >= 10)
  
  #combine all datasets into a single dataframe
  combined <- Input %>%
    left_join(WT_plasmid, by = c("site", "amplicon", "mutation", "mutation.type"),
              suffix = c(".pla", ".wt")) %>%
    left_join(P0, by = c("site", "amplicon", "mutation", "mutation.type")) %>% 
    rename("count.p0" = count) %>% 
    left_join(P1, by = c("site", "amplicon", "mutation", "mutation.type")) %>% 
    rename("count.p1" = count) %>% 
    left_join(P4, by = c("site", "amplicon", "mutation", "mutation.type")) %>% 
    rename("count.p4" = count) %>% 
    ungroup()
  
  #calculate frequencies (with a pseudocount added)
  #calculate frequencies based on amplicon
  #meaning using the equation: freqency_{i,amplicon} = \frac{read\:count_{i, amplicon}+1}{\displaystyle\sum\limits_{k \in amplicon}(read\:count_{k,amplicon}+1)}
  freq_amplicon <- combined %>% 
    group_by(amplicon) %>% 
    mutate(freq.pla = (count.pla+1)/(sum(count.pla) + length(count.pla)),
           freq.wt = (count.wt+1)/(sum(count.wt) + length(count.wt)),
           freq.p0 = (count.p0+1)/(sum(count.p0) + length(count.p0)),
           freq.p1 = (count.p1+1)/(sum(count.p1) + length(count.p1)),
           freq.p4 = (count.p4+1)/(sum(count.p4) + length(count.p4))) %>%
    #filter for frequency in input >= 6 * frequency in wt plasmid
    filter(freq.pla >= 6*freq.wt |
             mutation.type == "wildtype")
  
  enrichment_amplicon <- freq_amplicon %>%
    mutate(enrichment.p0 = (freq.p0)/(freq.pla),
           enrichment.p1 = (freq.p1)/(freq.pla),
           enrichment.p4 = (freq.p4)/(freq.pla))
  
  return(enrichment_amplicon)
}

############################## func 7 ###################################
# used in amino acid decrease section rather than amino acid preferences section
# function after "compute_enrichment", used to calculate fitness (log enrichment ratio)
# parallel to function "compute_nolog_fitness", which is used in entropy calculation
# compute fitness by replicate (instead of pre-, wt-, post-passage)
# compute fitness pla <- p0, pla <- p1, pla <- p4
compute_fitness_by_replicate <- function(enrichment_amplicon_file){
  
  #calculate mean enrichment of silent mutations
  amplicon_silent_mean <- enrichment_amplicon_file %>%
    filter(mutation.type == "silent") %>%
    group_by(amplicon) %>% 
    summarize(silent.mean.p0 = mean(log10(enrichment.p0)),
              silent.mean.p1 = mean(log10(enrichment.p1)),
              silent.mean.p4 = mean(log10(enrichment.p4)))
  
  fitness <- enrichment_amplicon_file %>% 
    left_join(amplicon_silent_mean, by = "amplicon") %>% 
    mutate(fitness.p0 = ifelse(mutation.type == "wildtype", log10(enrichment.p0),
                               log10(enrichment.p0) - silent.mean.p0),
           fitness.p1 = ifelse(mutation.type == "wildtype", log10(enrichment.p1),
                               log10(enrichment.p1) - silent.mean.p1),
           fitness.p4 = ifelse(mutation.type == "wildtype", log10(enrichment.p4),
                               log10(enrichment.p4) - silent.mean.p4))
  
  return(fitness)
}

############################## func 8 ###################################
# function that extracts the mutations that were present in the plasmid libraries
# but not P0
AA0count <- function(Plasmidfile, P0file){
  Mut_plasmid <- read_counts(Plasmidfile)
  WT_plasmid <- read_counts("WTPla")
  Mut_P0 <- read_counts(P0file)
  
  combined <- Mut_plasmid %>%
    left_join(WT_plasmid, by = c("site", "amplicon", "mutation", "mutation.type"),
              suffix = c(".mtpla", ".wtpla")) %>%
    left_join(Mut_P0, by = c("site", "amplicon", "mutation", "mutation.type")) %>%
    rename("count.mtp0" = count) %>%
    ungroup()
  
  freq_amplicon <- combined %>% 
    group_by(amplicon) %>% 
    mutate(freq.mtpla = (count.mtpla+1)/(sum(count.mtpla) + length(count.mtpla)),
           freq.wtpla = (count.wtpla+1)/(sum(count.wtpla) + length(count.wtpla)),
           freq.mtp0 = (count.mtp0+1)/(sum(count.mtp0) + length(count.mtp0))) %>%
    #filter for frequency in input >= 6 * frequency in wt plasmid
    filter(freq.mtpla >= 6*freq.wtpla |
             mutation.type == "wildtype")
  
  lethal <- filter(freq_amplicon, count.mtp0 == 0 & count.mtpla > 0)$mutation
  
  return(lethal)
}

############################## func 9 ###################################
# function after "read_counts", used to calculte the enrichment ratio for amino acid preference
compute_enrichment <- function(prepassage, WT_pla, postpassage){
  Mut_plasmid <- read_counts(prepassage)
  WT_plasmid <- read_counts(WT_pla)
  passaged <- read_counts(postpassage)
  
  #filter for counts >10 in plasmid libary
  Input <- Mut_plasmid %>%
    filter(count >= 10)
  
  #combine all datasets into a single dataframe
  combined <- Input %>%
    left_join(WT_plasmid, by = c("site", "amplicon", "mutation", "mutation.type"),
              suffix = c(".input", ".wt")) %>%
    left_join(passaged, by = c("site", "amplicon", "mutation", "mutation.type")) %>%
    rename("count.passaged" = count) %>%
    ungroup()
  
  #calculate frequencies (with a pseudocount added)
  #calculate frequencies based on amplicon
  #meaning using the equation: freqency_{i,amplicon} = \frac{read\:count_{i, amplicon}+1}{\displaystyle\sum\limits_{k \in amplicon}(read\:count_{k,amplicon}+1)}
  freq_amplicon <- combined %>% 
    group_by(amplicon) %>% 
    mutate(freq.input = (count.input+1)/(sum(count.input) + length(count.input)),
           freq.wt = (count.wt+1)/(sum(count.wt) + length(count.wt)),
           freq.passaged = (count.passaged+1)/(sum(count.passaged) + length(count.passaged))) %>%
    #filter for frequency in input >= 6 * frequency in wt plasmid
    filter(freq.input >= 6*freq.wt |
             mutation.type == "wildtype")
  
  enrichment_amplicon <- freq_amplicon %>%
    mutate(enrichment = (freq.passaged)/(freq.input))
  
  return(enrichment_amplicon)
}

############################## func 10 ###################################
# function after "compute_enrichment", used to calculate fitness (log enrichment ratio)
# parallel to function "compute_nolog_fitness", which is used in entropy calculation
compute_fitness <- function(enrichment_amplicon_file){
  
  #calculate mean enrichment of silent mutations
  amplicon_silent_mean <- enrichment_amplicon_file %>%
    filter(mutation.type == "silent") %>%
    group_by(amplicon) %>% 
    summarize(silent.mean = mean(log10(enrichment)))
  
  fitness <- enrichment_amplicon_file %>% 
    left_join(amplicon_silent_mean, by = "amplicon") %>% 
    mutate(fitness = ifelse(mutation.type == "wildtype", log10(enrichment),
                            log10(enrichment) - silent.mean))
  
  return(fitness)
}

############################## func 11 ###################################
# function after "compute_enrichment", nolog fitness is the same thing as enrichment ratio,
# but is used in entropy calculation
compute_nolog_fitness <- function(enrichment_amplicon_file){
  
  nolog_amplicon_silent_mean <- enrichment_amplicon_file %>%
    filter(mutation.type == "silent") %>%
    group_by(amplicon) %>%
    #no log10 here
    summarize(silent.mean = mean(enrichment))
  
  nolog_fitness <- enrichment_amplicon_file %>% 
    left_join(nolog_amplicon_silent_mean, by = "amplicon") %>% 
    #divide instead of substract because we got rid of log10
    mutate(fitness = ifelse(mutation.type == "wildtype", enrichment,
                            enrichment/silent.mean))
  
  return(nolog_fitness)
}

############################## func 12 ###################################
# function for entropy plotting
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

########################### helper df 3 #################################
# translate is almost the same as "codons", only difference is upper/lower case
translate <- tibble(
  codon = tolower(
    c("TTT","TTC",
      "TTA","TTG","CTT","CTC","CTA","CTG",
      "ATT","ATC","ATA",
      "ATG",
      "GTT","GTC","GTA","GTG",
      "TCT","TCC","TCA","TCG","AGT","AGC",
      "CCT","CCC","CCA","CCG",
      "ACT","ACC","ACA","ACG",
      "GCT","GCC","GCA","GCG",
      "TAT","TAC",
      "TAA","TAG","TGA",
      "CAT","CAC",
      "CAA","CAG",
      "AAT","AAC",
      "AAA","AAG",
      "GAT","GAC",
      "GAA","GAG",
      "TGT","TGC",
      "TGG",
      "CGT","CGC","CGA","CGG","AGA","AGG",
      "GGT","GGC","GGA","GGG")
  ),
  aa = c("Phe","Phe",
         "Leu","Leu","Leu","Leu","Leu","Leu",
         "Ile","Ile","Ile",
         "Met",
         "Val","Val","Val","Val",
         "Ser","Ser","Ser","Ser","Ser","Ser",
         "Pro","Pro","Pro","Pro",
         "Thr","Thr","Thr","Thr",
         "Ala","Ala","Ala","Ala",
         "Tyr","Tyr",
         "***","***","***",
         "His","His",
         "Gln","Gln",
         "Asn","Asn",
         "Lys","Lys",
         "Asp","Asp",
         "Glu","Glu",
         "Cys","Cys",
         "Trp",
         "Arg","Arg","Arg","Arg","Arg","Arg",
         "Gly","Gly","Gly","Gly")
)

############################## func 13 ###################################
# readin fasta file: a fasta file set with multiple sequences
# output have one nucleotide/amino acid in every column, the number of column equals to the number of
# residues (758); each row is a sample, rownames are the Isolate_Id in the fasta file.
readMSA <- function(file_path){
  t <- readLines(file_path)
  id_index <- which(grepl(">", t))
  ids <- gsub(">", "", t[id_index])
  nseqs <- length(id_index)
  
  seqs <- numeric(nseqs)
  for(i in 1:(nseqs-1)){
    seqs[i] <- paste(t[(id_index[i]+1):(id_index[i+1]-1)], collapse = "")
  }
  seqs[length(id_index)] <- paste(t[(id_index[nseqs]+1):length(t)], collapse = "")
  
  seq_mat <- matrix(unlist(strsplit(seqs, "")), ncol = nchar(seqs[1]), byrow = TRUE)
  rownames(seq_mat) <- ids
  
  return(seq_mat)
}

############################## func 14 ###################################
# get weighted entropy: this function works for one site, meaning the input 'mat' is a vector
# @mat: should be each column/residue of 'aa_mat_pre09' or 'aa_mat_post09'
# @df: the dataframe specifies the year of each sequence, 'df_pre09' or 'df_post09'
get_ent_weighted <- function(mat, df, weights=NULL){
  groups <- sort(unique(df$group))
  
  group_freqs <- list()
  for(i in 1:length(groups)){
    mat_i <- mat[df[df$group == groups[i],]$Isolate_Id,]
    group_freqs[[i]] <- as.list(prop.table(table(mat_i)))
  }
  
  group_freqs_bind <- data.table::rbindlist(group_freqs, fill = TRUE)
  group_freqs_bind[is.na(group_freqs_bind)] <- 0
  
  if(is.null(dim(weights))){
    mean_freqs <- colMeans(group_freqs_bind)
  }
  else{
    mean_freqs <- colSums(group_freqs_bind*weights$weight)
  }
  
  return(sum(-log(mean_freqs)*mean_freqs))
}

############################## func 15 ###################################
# pair with func 'get_ent_weighted', to iterate over all 758 sites
# @mat: aa_mat_pre09, aa_mat_post09
# @df: df_pre09, df_post09
get_ent_all_sites <- function(mat, df){
  shannon_diversity <- vector(length = ncol(mat))
  for(i in 1:ncol(mat)){
    matrix_i <- as.matrix(mat[,i])
    rownames(matrix_i) <- rownames(mat)
    shannon_diversity[i] = get_ent_weighted(matrix_i, df)
  }
  
  return(shannon_diversity)
}

############################## func 16 ###################################
# Calculate the minimum nucleotide changes needed for an amino acid to become another one
# setup codon distance table (need 'codons' from main)
# str1 and str2 are vectors that can contain more than one string
# str1 and str2 are codons
compare_str <- function(str1, str2){
  comparison_set <- data.frame()
  
  for (each_piece1 in seq(1, length(str1))){
    for (each_piece2 in seq(1, length(str2))){
      string1 <- str_split(str1, "")[[each_piece1]]
      string2 <- str_split(str2, "")[[each_piece2]]
      codon_diff <- sum(string1 != string2)
      
      vr <- c(str1[[each_piece1]], str2[[each_piece2]], codon_diff)
      comparison_set <- rbind(comparison_set, vr)
    }
  }
  colnames(comparison_set) <- c("Str1", "Str2", "CodonDiff")
  
  return(comparison_set)
}

############################## func 17 ###################################
# to modify PDB file, using entropy data in replacement of b-factors
# @filename: quoted, pdb file name
# @entropy: a df with two (or more) columns; colnames(resno, avg.entropy)
# This is only an example function. The PDBs were in fact modified by executing the code line by line,
# because the numbering of chains and resnos for each PDB file is different.
modify_PDB <- function(filename, entropy){
  polymerase <- read.pdb(paste0(indir, "pdb/", filename), maxlines = -1, multi = FALSE, rm.insert = FALSE,
                         rm.alt = FALSE, ATOM.only = FALSE, hex = FALSE, verbose = TRUE)
  
  replace <- entropy %>% 
    mutate(chain = "B") %>% 
    select(chain, resno, average.entropy)
  
  pb1_modified <- merge(polymerase$atom, replace, by = c("chain", "resno"), all.x = TRUE)
  pb1_modified <- pb1_modified %>% 
    mutate(b = case_when(!is.na(average.entropy) ~ formatC(average.entropy, digits = 2, format = "f"), 
                         T ~ formatC(0, digits = 2, format = "f"))) %>% 
    select(-average.entropy)
  
  pb1_modified <- pb1_modified[order(pb1_modified$eleno), ]
  
  polymerase$atom <- pb1_modified
  
  #check
  # pb1 <- atom.select(polymerase, chain="B")
  # pb1_pdb <- trim.pdb(polymerase, pb1)
  # head(pb1_pdb$atom)
  
  write.pdb(polymerase, file=paste0(indir, "pdb/modified_", filename))
  #manually added terminal lines between different chains; manually added the beginning text.
}

############################## func 18 ###################################
# to get dataframe with rescaled enrichment ratio/fitness and specific sites,
# then transform the dataframe to matrix for ggseqlogo
# @sites: numbers can be c(1, 2, 3) or 1:3 or number vectors
# @entropy_calculation_df: the intermediate dataframe before 'entropy': entropy_calculation, it
# has site, mutation, avg.nolog.fitness, nolog.fitness.rescaled, and entropy. We need 'site' and
# 'nolog.fitness.rescaled'
get_logo_mat <- function(sites, entropy_calculation_df){
  logo <- entropy_calculation_df %>% 
    filter(site %in% sites) %>% 
    select(site, mutation, avg.nolog.fitness, nolog.fitness.rescaled, entropy) %>% 
    arrange(site) %>% 
    mutate(mutation.wt = substr(gsub('\\d','', mutation), 1, 3),
           mutation.mt = substr(gsub('\\d','', mutation), 4, 6))
  
  logo$wt <- aa321(toupper(logo$mutation.wt))
  logo$mt <- aa321(toupper(logo$mutation.mt))
  # logo$mt[logo$mt=="X"]<-"*"
  
  logo_mat <- spread(select(logo, site, mt, nolog.fitness.rescaled), 
                     key = as.numeric(site), value = nolog.fitness.rescaled)
  
  logo_mat <- replace(logo_mat, is.na(logo_mat), 0)
  
  logomatrix <- as.matrix(logo_mat[,-1])
  rownames(logomatrix)<-logo_mat$mt
  
  return(logomatrix)
}