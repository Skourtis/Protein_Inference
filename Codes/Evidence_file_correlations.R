#ProteomeHD exploration#
#Summarise abundance for PGs
#Do peptides that anticorrelate have specific characteristics, could the be put into a "black"list of quantification, do they have more modfications?
columns_to_read <- c("Sequence", "Proteins", "Experiment", "Protein group IDs","Leading Razor Protein",
                     "Ratio H/L normalized", "Ratio M/L normalized", "Ratio H/M normalized",
                     "Reverse", "Potential contaminant")
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
interesting_protein <- "P78371"
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")
# SQL_loading <-  read.csv.sql(here::here("Datasets","Raw","ProteomeHD", "evidence.txt"), sql = "select * from file where Proteins == 'P55011;P55011-3'")
testing <- readr::read_delim_chunked(file_input, delim = "\t",
                        callback = DataFrameCallback$new(function(x, pos) select(x,all_of(evidence_selection )) %>% 
                                                             janitor::clean_names() %>% 
                                                             subset(leading_razor_protein == interesting_protein))) %>% 
    setDT
sampling_evidence <- readr::read_delim_chunked(file_input, delim = "\t",chunk_size = 100000,
                                     callback = DataFrameCallback$new(function(x, pos) select(x,all_of(c("Proteins",#"Leading Proteins" ,
                                                                                                         "Leading Razor Protein"))) %>% 
                                                                          distinct())) %>% distinct() %>% 
    janitor::clean_names() %>% setDT
# inte sampling_evidence %>% subset()
# VRoom_loading <- vroom::vroom(here::here("Datasets","Raw","ProteomeHD", "evidence.txt"), col_select  = evidence_selection)
# interesting_protein_df <- VRoom_loading %>% 
#     subset(Proteins == interesting_protein)
Protein_groups <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = c("Intensity", "Protein IDs","Majority protein IDs")) %>%  
     janitor::clean_names() %>% 
    mutate(intensity = as.numeric(intensity)) 
# %>% 
    arrange(-intensity) %>% 
    head(100)
 testing_identical_peptides <- fread(input = here::here("Datasets","Raw","ProteomeHD", "evidence.txt"), nrows =  1000) %>% 
     janitor::clean_names()
# testing_identical_peptides %>% subset(experiment == "g4_PX441_F5" & sequence == "ADLINNLGTIAK") %>% View()
most_peptides <-testing_identical_peptides[, .N, by=.(leading_razor_protein)][order(-rank(N)),][1,1] %>% deframe()
most_peptides_evidence_leading <- testing_identical_peptides[leading_razor_protein == most_peptides,]
 most_peptides_evidence <- most_peptides_evidence_leading[, .(median_intensity = median(ratio_m_l_normalized,na.rm = T)),c("sequence","experiment")]


    most_peptides_evidence <- testing[, .(median_intensity = median(as.numeric(ratio_m_l_normalized),na.rm = T)),c("sequence","experiment")]
testing2 <- most_peptides_evidence[!is.nan(median_intensity),c("sequence","median_intensity","experiment")] %>% 
    dcast(.,
          experiment~sequence ,
          value.var = c("median_intensity"))
selected_peptides <- colSums(is.na(testing2)) < (nrow(testing2)*0.5 )
testing2 <- testing2[,..selected_peptides]
combos <- combinat::combn(testing2[,!("experiment")] %>% colnames(),2,simplify = F) %>% 
    set_names(.,purrr::map_chr(.x = .,~paste(.x[1],.x[2], sep = "_")))
unique_peptides <- testing %>% 
    subset(proteins == leading_razor_protein) %>% pull(sequence) %>% unique()
Peptide_correlation <- function(peptides_pair){
    # peptides_pair = combos[[1]]
    if( (testing2[,..peptides_pair] %>% 
        na.omit() %>% nrow())>30){
         # cor(testing2 %>% pull(peptides_pair[1]),testing2 %>% pull(peptides_pair[2]),method = "pearson", use = "pairwise")
        WGCNA::bicor(testing2 %>% pull(peptides_pair[1]),testing2 %>% pull(peptides_pair[2]),use = "pairwise.complete.obs")
    }else{
            NaN
        }
}
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
