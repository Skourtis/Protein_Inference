pacman::p_load(tidyverse,data.table)
#Recreating MQ Protein Intensities
#ProteomeHD exploration#
#Summarise abundance for PGs
#Do peptides that anticorrelate have specific characteristics, could the be put into a "black"list of quantification, do they have more modfications?
#Are all peptide evidence used for quantificaiton?
# MULTI-MSMS:     label pair detected
# - MULTI-SECPEP  identification through second peptide algorithm
# - MULTI-MATCH    identification based on match between runs
# - ISO-MSMS:         isotope cluster detected but no label pair
# - MSMS:                no MS1 label state detected

# https://groups.google.com/g/maxquant-list/c/06JU00WWfR8/m/SsDGgaWICQAJ
# it depends on the TYPE in evidence.
# But it is very complicate.
# If there are many rows with multiMSMS, than only these are used for Protein ratio.
# This will change the Protein ratio H/L in direction to 1, because if the ratio of H/L is far from 1 the probability that there is only an identified MSMS of one labeling state is high and these 
# peptide ratio are not taken into account for the protein ratio. This explains many of my observations.
# 
# If there are no or only a small number of rows with multiMSMS than ISO-MSMS rows are also included and the Ratio H/L iso-count  is >0 in ProteinGroups.
# I never found a ratio H/L for TYPE= MSMS.
# 
# 1. Which ratios are used for calculation of Protein ratio H/L ? Can the user manipulate which peptide ratio are used ?
#     
#     But there comes up the next questions:
#     2. Why some rows are TYPE= MSMS and other TYPE= ISO-MSMS. Looking at the same Peptide Sequence - It can not be explained by the PEP value.
# Looking in raw-data in both cases is there an Isotopic pattern with the corresponding charge state recognized by QualBrowser.
# 
# 3. Why are there many rows for the same unmodified peptide from the same raw-file ?
#     It seems me that each identified MSMS produces one row in Evidence table.
# That means additional rows are produced by :
#     - different charge states
# - repeated measurements of one charge state (Peak tailing longer than exclusion time)
# - during measurement wrong recognized monoisotopic masses

###PROBLEMS

##Proteins which have NaN in HL and values in ML, but do not use the values in ML
# Q969P0-3;Q969P0 in g3_px597_b3_spr

#But also the opposite happens where proteins have NaN in evidence HL, and they use the ML value in PG
# A0A183 in g1_px1406_gm18858 experiment


columns_to_read <- c("Sequence", "Proteins", "Experiment",# "Protein group IDs",
                     "Leading Razor Protein","id",
                     "Ratio H/L normalized", "Ratio M/L normalized",# "Ratio H/M normalized" ,"Peptide ID",
                     # "Reverse", "Potential contaminant",
                     "Type")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")

# Read in data
DT <- fread(file_input, select = columns_to_read)

#it seems that for mot cases, if the HL is NA and the ML is there,
#then the ML is the real HL and when used agrees with the PG ratio
#This is especially true for proteins with a single evidence
DT[,`Ratio H/L normalized`:= fifelse(is.na(`Ratio H/L normalized`)  #& (!is.nan(`Ratio H/L normalized`))
                                     ,`Ratio M/L normalized`,`Ratio H/L normalized`)]

## DT[,`Ratio H/L normalized`:= fifelse(is.nan(`Ratio M/L normalized`),NaN,`Ratio H/L normalized`)]
DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA

DT_cols <- fread(file_input, nrows = 1)

#loading the proteins to get the evidence for a protein
Protein_groups <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = c("Intensity", "Protein IDs","Majority protein IDs", "Evidence IDs", "Peptide IDs"))

`%out%` <- negate(`%in%`)
Protein_evidence <- Protein_groups %>%
  janitor::clean_names() %>% 
  dplyr::select(c(majority_protein_i_ds,evidence_i_ds)) %>% 
  separate_rows(evidence_i_ds) %>% 
  rename(id = evidence_i_ds) %>% 
  setDT() 
Protein_evidence[,id := as.numeric(id)]

Protein_groups_empty <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), nrows = 2)
cols_of_interest <- colnames(Protein_groups_empty) %>% 
  str_subset("Majority|Ratio H/L normalized ")
cols_Type <- colnames(Protein_groups_empty) %>% 
  str_subset("Majority|Ratio H/L normalized |H/L type| H/L iso-count ")

#loading the proteins to get both the ratios but also the Type and how many IsoMSMS were used for quantification

Long_protein_groups <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = cols_of_interest)%>% 
  janitor::clean_names() %>%  
  dplyr::select(matches("majority|ratio_h_l_normalized_")) %>%
  pivot_longer(cols  = contains("ratio"),names_to = "Experiment", values_to = "Ratio" ) %>% 
  subset(!is.na(Ratio) )%>% 
  setDT()
Long_protein_groups <- Long_protein_groups[,Experiment:= str_remove_all(Experiment,"ratio_h_l_normalized_")]
Long_protein_groups_type <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = cols_Type)%>% 
  janitor::clean_names() %>%  
  dplyr::select(matches("majority|ratio_h_l_normalized_|h_l_type_|iso_count")) %>% 
  setDT()
Ratio <- Long_protein_groups_type %>% 
  dplyr::select(matches("majority|ratio_h_l_normalized_")) %>% 
  melt(id.vars = c("majority_protein_i_ds"),
       measure.vars = patterns("^ratio_h_l_normalized"),
       variable.name = "Experiment", value.name = "h_l_ratio") %>% 
  .[,Experiment := str_remove_all(Experiment,"ratio_h_l_normalized_") ]
Type <- Long_protein_groups_type %>% 
  dplyr::select(matches("majority|ratio_h_l_type_")) %>% 
  melt(id.vars = c("majority_protein_i_ds"),
       measure.vars = patterns("^ratio_h_l_type"),
       variable.name = "Experiment", value.name = "type") %>% .[,Experiment := str_remove_all(Experiment,"ratio_h_l_type_") ]
Iso <- Long_protein_groups_type %>% 
  dplyr::select(matches("majority|iso_count")) %>% 
  melt(id.vars = c("majority_protein_i_ds"),
       measure.vars = patterns("^ratio_h_l_iso"),
       variable.name = "Experiment", value.name = "iso_count") %>% 
  .[,Experiment := str_remove_all(Experiment,"ratio_h_l_iso_count_") ]

# Counting number of ISOMSMS used for a protein
Long_protein_groups_type <- Ratio[Type, on = .(majority_protein_i_ds==majority_protein_i_ds, Experiment== Experiment)
][Iso, on = .(majority_protein_i_ds==majority_protein_i_ds, Experiment== Experiment)
][,rm_iso := fifelse(iso_count == 0,T,F)]

testing_evidence_sum <- DT[Protein_evidence, on = .(id = id)
][,leading_majority := str_remove_all(majority_protein_i_ds,";[:graph:]*$")
][`Leading Razor Protein` ==  leading_majority,
  ][is.finite(`Ratio H/L normalized`)
]

#cleaning the names of experiments
fixed_exp_names <- data.table(Experiment_old = testing_evidence_sum$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]
testing_evidence_sum <- testing_evidence_sum[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]
testing_evidence_sum <- setnames(testing_evidence_sum, c("i.Experiment"), c("Experiment"))

#removing ISOMSMS not used for quantification and counting valid evidence per protein
testing_evidence_sum <- testing_evidence_sum[Long_protein_groups_type[is.finite(h_l_ratio),], on = .(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)]
testing_evidence_sum <- testing_evidence_sum[!(rm_iso == T & Type == "ISO-MSMS")][!is.na(`Ratio H/L normalized`)]
evidence_n <-testing_evidence_sum[,.N,by = .(Experiment,majority_protein_i_ds)] 

#median of valid evidence
testing_evidence_sum<-testing_evidence_sum[,.(Prot_HL=median(`Ratio H/L normalized`,na.rm = T)), by = .(majority_protein_i_ds,Experiment)
]

#combine with PG ratio
PG_summary_comparison <- testing_evidence_sum[ Long_protein_groups, on = .(Experiment = Experiment,majority_protein_i_ds = majority_protein_i_ds)][
  , `:=`(Prot_HL = round(Prot_HL,3),
         Ratio = round(Ratio,3))
]

#removing proteins which didn't arise from median (advanced ratio calculation)
PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  
][evidence_n,on = .(majority_protein_i_ds = majority_protein_i_ds,Experiment= Experiment)
][,`:=`(difference = Prot_HL-Ratio)
][Long_protein_groups_type, on = .(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)
][type == "Median",]#[is.na(Prot_HL) |is.na(Ratio),]


#plots
PG_summary_comparison[Identical == T,N] %>% .[.<40] %>% 
  hist(breaks = 1000, main = "Hist of evidence number of Correct PG ratio", xlim = c(0,40))
PG_summary_comparison[Identical == F,N] %>% .[.<40] %>% hist(breaks = 1000, main = "Hist of evidence number of Incorrect PG ratio", xlim = c(0,40))

PG_summary_comparison[!(is.na(Prot_HL)|is.na(Ratio))] %>%  
  ggplot(aes(x = log2(Prot_HL), y = log2(Ratio), colour= fifelse(N%%2 == 0,"Even","Odd")))+
  geom_point(alpha = 0.3)+labs(title="Recalculating ProteinGroup Silac Ratio from Evidence file",
                               x ="Ratio Recalculated from Evidence file", y = "Ratio from ProteinGroup",colour = "Number of Evidence for PG in Experiment")
  
ggsave(here::here("Output","Silac_Ratio_PG_reconstruction.png"))


experiment = "g1_px1194_pca1_1"
protein_id = "Q9Y3C5"
DT[,Experiment := tolower(Experiment)]
 DT[str_detect(Proteins,protein_id) &  str_detect(Experiment,experiment)# & Type == "MULTI-MSMS"
    ,] 
# testing_evidence_sum[str_detect(majority_protein_i_ds,protein_id) &  str_detect(Experiment,experiment)]
# evidence_ids <- DT[str_detect(Proteins,protein_id) &  str_detect(Experiment,experiment),id]
# Protein_evidence[id %chin% evidence_ids]



