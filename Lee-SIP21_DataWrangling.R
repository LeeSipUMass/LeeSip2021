### setup: install/load packages, set WD

install.packages("tidyverse")
library(tidyr)
library(dplyr)       
library(tidyverse)
setwd("R Materials/")

### import orthologs
# fetched orthologs through browser (not R), downloaded to WD
# import homologues gene set [Ensembl BioMart]
PFCEP_orthologs <- read_tsv("PFCEP_Orthologs_withParalogs_3species.txt")

# rename columns for ease of coding, then (pipe) remove duplicate rows, then remove rows with NAs
PFCEP_orthologs <- rename(PFCEP_orthologs, 
       Hum = "Human_(GRCh38.p13)_Gene_stable_ID", # human gene ID
       Chi = "Chimpanzee_(Pan_tro_3.0)_gene_stable_ID", 
       Gor = "Gorilla_(gorGor4)_gene_stable_ID"
        ) %>% distinct() %>% drop_na()
dim(PFCEP_orthologs)

# keep only 1:1 paralogs, then remove paralog columns
PFCEP_orthologs <- filter(PFCEP_orthologs, Chi_paralogs == "ortholog_one2one") %>% select (-Chi_paralogs)
PFCEP_orthologs <- filter(PFCEP_orthologs, Gor_paralogs == "ortholog_one2one") %>% select (-Gor_paralogs)

dim(PFCEP_orthologs)

# create species pairs of orthologs with species genes ("Hum") in first column
Hum_Hum <- select(PFCEP_orthologs, Hum)
Hum_Chi <- select(PFCEP_orthologs, Hum, Chi) 
Hum_Gor <- select(PFCEP_orthologs, Hum, Gor)

### import sample count data, then rename columns for ease of joining
# txt files have two columns:
  # gene: species-specific gene ID
  # count: data from RNA-Seq experiments
Hum1 <- read_tsv("htseq_Hum1.txt") %>% rename(Hum = "gene", Hum1 = "count") 
Hum2 <- read_tsv("htseq_Hum2.txt") %>% rename(Hum = "gene", Hum2 = "count")
Hum3 <- read_tsv("htseq_Hum3.txt") %>% rename(Hum = "gene", Hum3 = "count")
Chi1 <- read_tsv("htseq_Chi1.txt") %>% rename(Chi = "gene", Chi1 = "count") 
Chi2 <- read_tsv("htseq_Chi2.txt") %>% rename(Chi = "gene", Chi2 = "count")
Chi3 <- read_tsv("htseq_Chi3.txt") %>% rename(Chi = "gene", Chi3 = "count")
Gor1 <- read_tsv("htseq_Gor1.txt") %>% rename(Gor = "gene", Gor1 = "count") 

# match sample count data to orthologous human genes, then select out NHP gene ID
Hum1 <- inner_join(Hum_Hum, Hum1, by = "Hum")
Hum2 <- inner_join(Hum_Hum, Hum2, by = "Hum")
Hum3 <- inner_join(Hum_Hum, Hum3, by = "Hum")
Chi1 <- inner_join(Hum_Chi, Chi1, by = "Chi") %>% select (-Chi)
Chi2 <- inner_join(Hum_Chi, Chi2, by = "Chi") %>% select (-Chi)
Chi3 <- inner_join(Hum_Chi, Chi3, by = "Chi") %>% select (-Chi)
Gor1 <- inner_join(Hum_Gor, Gor1, by = "Gor") %>% select (-Gor)

### combine counts within species
Hum_counts <- inner_join(Hum1, Hum2) %>% inner_join(Hum3)
Chi_counts <- inner_join(Chi1, Chi2) %>% inner_join(Chi3) 
Gor_counts <- Gor1

### combine counts across species
PFCEP_counts <- Hum_counts
PFCEP_counts <- inner_join(PFCEP_counts, Chi_counts, by = "Hum")
PFCEP_counts <- inner_join(PFCEP_counts, Gor_counts, by = "Hum")

dim(PFCEP_counts)
# View(PFCEP_counts)

# rename gene column
PFCEP_counts <- rename(PFCEP_counts, gene = "Hum")

# export final dataset
write_tsv(PFCEP_counts, file = "PFCEP_counts_3species.txt")