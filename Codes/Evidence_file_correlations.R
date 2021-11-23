#ProteomeHD exploration#
#Summarise abundance for PGs
#Do peptides that anticorrelate have specific characteristics, could the be put into a "black"list of quantification, do they have more modfications?
columns_to_read <- c("Sequence", "Proteins", "Experiment", "Protein group IDs","Leading Razor Protein",
                     "Ratio H/L normalized", "Ratio M/L normalized", "Ratio H/M normalized",
                     "Reverse", "Potential contaminant")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")

# Read in data
DT <- fread(file_input, select = columns_to_read)
# Drop reverse and contaminant peptides
DT <- DT[ Reverse != "+" & `Potential contaminant` != "+" ]

# What's going with the H/L, M/L, H/M ratios? W
DT[, lapply(.SD, function(x){ sum(!is.na(x)) }), .SDcols = c("Ratio H/L normalized", "Ratio M/L normalized", "Ratio H/M normalized") ]
# -> 16 million M/L values when most experiments didn't even have a medium SILAC label?

# Perhaphs
DT[, lapply(.SD, function(x){ sum(!is.na(x)) }), 
   .SDcols = c("Ratio H/L normalized", "Ratio M/L normalized",  "Ratio H/M normalized"), 
   by = gsub("_.+", "", Experiment)]
# -> This shows H/L and H/M values are only assigned for parameter group g3, which is the triple labelled group
# MQ must internally re-name these 

# For now, just work with the peptides that have an M/L ratio
DT <- DT[ !is.na(`Ratio M/L normalized`) ]

# Drop unnecessary columns
DT[, c("Ratio H/L normalized", "Ratio H/M normalized", "Reverse", "Potential contaminant") := NULL ]

# How many different peptides & experiments are there, and how often have they be seen?
DT[, .N, Sequence][ order(-N) ]
DT[, .N, Experiment][ order(-N) ]

# What about proteins?
DT[, .N, Proteins][ order(-N) ]             # ~22k proteins?
DT[, .N, `Protein group IDs`][ order(-N) ]  # ~15k final protein groups? (this must be before some QC cut-off!)

# Note that some protein groups contain multiple proteins
DT[ `Protein group IDs` == "8265", .N, .(Proteins, `Protein group IDs`)]

# And some proteins are assigned to multiple proteinGroups
DT[ Proteins == "Q5U651", .N, .(Proteins, `Protein group IDs`)]

# How many observations per peptide & experiment (per proteinGroup or per protein? I start with Proteins)? 
DT[, .N, .(Sequence, Proteins, Experiment)][ order(-N) ]

# For each peptide, get median M/L ratio per Experiment
med_DT <- DT[, .(median_ratio = median(`Ratio M/L normalized`, na.rm = TRUE)), .(Sequence, Proteins, Experiment)]

# Exlude peptides that were seen in less than 25 experiments, since we can't create a solid correlation profile
exlude_peptides <- med_DT[, .N, Sequence][ N <= 25, Sequence ]
med_DT <- med_DT[ !Sequence %in% exlude_peptides]

# Calculate correlations among peptides for each protein
pb <- txtProgressBar(min = 0, max =length(unique(med_DT$Proteins)), style = 3)
correlations <- data.table()
#maybe instead of cycling through Proteins, we should cycle through Leading Razor proteins,
#Because the Proteins and LEading proteins are treated differently here but will be used for the quant of the same Protein
for(i in 1:length(unique(med_DT$Proteins))){
    tmp_protein <- unique(med_DT$Proteins)[i]
    tmp <- med_DT[ Proteins == tmp_protein ]
    # Should we log 2 the median Ratio?  now the positive FC might have a larger influence on corr than negative
    # tmp[, median_ratio := log2(median_ratio)]
     # tmp %>% ggplot(aes(x  = Experiment, y = median_ratio, group = Sequence, colour = Sequence))+geom_point()+ggtitle(tmp_protein)
    tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
    tmp$Experiment <- NULL
    #finding and removing correlations with less than 15 common datapoints
    missingness_tmp <- crossprod(is.na(tmp))>15
    #using suggested bicor
    tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")
    tmp[missingness_tmp == FALSE] <- NaN
    tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))
    tmp <- tmp[, .( Peptide1 = as.character(Var1), Peptide2 = as.character(Var2), PCC = value ) ]
    tmp <- tmp[ Peptide1 > Peptide2 ][!is.nan(PCC)]
    #if using leading razor, then we can then look back at a simplified evidence table to see if the peptides were unique or shared in Proteins
    tmp[, Protein := tmp_protein ]
    correlations <- rbind( correlations, tmp)
    setTxtProgressBar(pb, i)
}

evidence_selection <- c("Sequence",
                        # "Length",
                        # "Modifications","Missed cleavages",
                        "Proteins",#"Leading Proteins" ,
                        "Leading Razor Protein",
                        # "Potential contaminant","Reverse","Intensity","Intensity L", "Intensity M","Intensity H",
                        "Ratio M/L normalized",
                        # "Ratio H/L normalized","Ratio H/M normalized", "PEP",
                        "Experiment")
sampling_evidence <- readr::read_delim_chunked(file_input, delim = "\t",chunk_size = 100000,
                                     callback = DataFrameCallback$new(function(x, pos) select(x,all_of(c("Proteins",#"Leading Proteins" ,
                                                                                                         "Leading Razor Protein"))) %>% 
                                                                          distinct())) %>% distinct() %>% 
    janitor::clean_names() %>% setDT
Protein_groups <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = c("Intensity", "Protein IDs","Majority protein IDs")) %>%  
     janitor::clean_names() %>% 
    mutate(intensity = as.numeric(intensity)) 
Protein_summary <- purrr::map_dbl(combos,Peptide_correlation) %>% enframe("Pair","Correlation") %>% 
    separate(Pair, into = c("peptide1", "peptide2")) %>% 
    mutate(Pair_type = case_when(
        peptide1 %in% unique_peptides &  peptide2 %in% unique_peptides ~ "Pair_of_uniques",
        peptide1 %in% unique_peptides & !(peptide2 %in% unique_peptides) ~ "One_unique",
       !(peptide1 %in% unique_peptides) & peptide2 %in% unique_peptides ~ "One_unique",
       !(peptide1 %in% unique_peptides) & !(peptide2 %in% unique_peptides) ~ "None_unique",
       TRUE ~"Other"),
       leading_razor_protein = interesting_protein)
Peptide_abundances <- melt(testing2, id.vars = "experiment",
             measure.vars = colnames(testing2[,!("experiment")]),
     variable.name = "peptide"  )
Peptide_abundances[, median(value,na.rm = T), by = .(experiment)]
Peptide_abundances <- Peptide_abundances[Peptide_abundances[, median(value,na.rm = T), by = .(experiment)], on = .(experiment = experiment)]
list(Peptide_abundances)

Peptide_abundances %>% 
    ggplot(aes(x = reorder(experiment,V1), y = value, group = peptide, colour  = peptide))+
    geom_point()+
    geom_line()
Protein_summary %>% ggplot(aes(x = Pair_type, y = Correlation ))+
    geom_boxplot()+
    geom_point()+
    ggtitle(glue::glue(interesting_protein, " Peptide_Correlation"))
