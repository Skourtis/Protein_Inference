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



columns_to_read <- c("Sequence", "Proteins", "Experiment",# "Protein group IDs",
                     "Leading Razor Protein","id",
                     "Ratio H/L normalized", "Ratio M/L normalized",# "Ratio H/M normalized" ,"Peptide ID",
                     # "Reverse", "Potential contaminant",
                     "Type")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")

# Read in data
DT <- fread(file_input, select = columns_to_read)
DT[,`Ratio H/L normalized`:= if_else(is.na(`Ratio H/L normalized`) & !is.nan(`Ratio H/L normalized`),`Ratio M/L normalized`,`Ratio H/L normalized`)]
DT[,`Ratio H/L normalized`:= if_else(is.nan(`Ratio M/L normalized`),NaN,`Ratio H/L normalized`)]

DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA
DT_cols <- fread(file_input, nrows = 1)

# Drop reverse and contaminant peptides
#DT <- DT[ Reverse != "+" & `Potential contaminant` != "+" ]

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
med_DT <- DT[, .(median_ratio = median(`Ratio M/L normalized`, na.rm = TRUE)), .(Sequence, Proteins, Experiment, `Leading Razor Protein`)]
pep_PG <- unique(med_DT[,.(Sequence, Proteins, `Leading Razor Protein`)])
pep_PG <- pep_PG[,unique:= if_else( Proteins ==`Leading Razor Protein`,"unique","non_unique" )][unique == "unique",.(Sequence)][["Sequence"]]

# Exlude peptides that were seen in less than 25 experiments, since we can't create a solid correlation profile
exlude_peptides <- med_DT[, .N, Sequence][ N <= 25, Sequence ]
med_DT <- med_DT[ !Sequence %in% exlude_peptides]

# Calculate correlations among peptides for each protein
pb <- txtProgressBar(min = 0, max =length(unique(med_DT$`Leading Razor Protein`)), style = 3)
correlations <- data.table()
#maybe instead of cycling through Proteins, we should cycle through Leading Razor proteins,
#Because the Proteins and LEading proteins are treated differently here but will be used for the quant of the same Protein
for(i in 1:length(unique(med_DT$`Leading Razor Protein`))){
    tmp_protein <- unique(med_DT$`Leading Razor Protein`)[i]
    tmp <- med_DT[ `Leading Razor Protein` == tmp_protein ]
    # Should we log 2 the median Ratio?  now the positive FC might have a larger influence on corr than negative
    tmp[, median_ratio := log2(median_ratio)]
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
pep1_unique <- correlations$Peptide1 %in% pep_PG + correlations$Peptide2 %in% pep_PG
correlations[,pair_type := fcase(pep1_unique == 2, "both_unique",
                                 pep1_unique ==1, "one_unique",
                                 pep1_unique ==0 , "none_unique")]

intersect(correlations[pair_type  == "one_unique" ,"Protein"][["Protein"]] %>% unique(),
          correlations[pair_type  == "both_unique" ,"Protein"][["Protein"]] %>% unique())
correlations[,.N,by = Protein] %>% arrange(-N)
matrix_1 <-  dcast(correlations[Protein  == "Q14204",], Peptide1~Peptide2, value.var = "PCC") %>% column_to_rownames("Peptide1") 
matrix_1 <- matrix_1[colnames(matrix_1),colnames(matrix_1)] 
matrix_1[is.na(matrix_1)] <- 0
matrix_1 = matrix_1+t(matrix_1)
pheatmap::pheatmap(matrix_1, cluster_rows = T,cluster_cols = T)
df.dist <- as.dist(matrix_1)
correlations[Protein  == "Q15149",] %>% 
    ggplot(aes(x = pair_type, y = PCC))+
    geom_boxplot()+
    geom_point()
Protein_groups_empty <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), nrows = 2)
Protein_groups <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = c("Intensity", "Protein IDs","Majority protein IDs", "Evidence IDs", "Peptide IDs"))
                        
# evidence_selection <- c("Sequence",
#                         # "Length",
#                         # "Modifications","Missed cleavages",
#                         "Proteins",#"Leading Proteins" ,
#                         "Leading Razor Protein",
#                         # "Potential contaminant","Reverse","Intensity","Intensity L", "Intensity M","Intensity H",
#                         "Ratio M/L normalized",
#                         # "Ratio H/L normalized","Ratio H/M normalized", "PEP",
#                         "Experiment")
# sampling_evidence <- readr::read_delim_chunked(file_input, delim = "\t",chunk_size = 100000,
#                                      callback = DataFrameCallback$new(function(x, pos) select(x,all_of(c("id" ,"Peptide ID",
#                                                                                                          columns_to_read))) %>% 
#                                                                           distinct())) %>% distinct() %>% 
#     # janitor::clean_names() %>% 
#     setDT
`%out%` <- negate(`%in%`)
Protein_evidence <- Protein_groups %>%
    janitor::clean_names() %>% 
    dplyr::select(c(majority_protein_i_ds,evidence_i_ds)) %>% 
    separate_rows(evidence_i_ds) %>% 
    rename(id = evidence_i_ds) %>% 
    setDT() 
Protein_evidence[,id := as.numeric(id)]

evidence_used <-  Protein_evidence$id %>% unique()

#apparently all evidence entries from the evidence file are being used
#even if contaminant or reverse
evidence_not_used <- DT[id %out% evidence_used,]

#this is explained as the reverse and cont evidence map to uniprots with CON_ or REV_
reverse_cont_evidence <-  DT[ Reverse == "+", ]
Protein_evidence[id == 20175,]
# DT[is.nan(DT)] <- NA
# DT <- DT[,`Ratio H/L normalized`:= if_else(is.na(`Ratio H/L normalized`),`Ratio M/L normalized`,`Ratio H/L normalized`)]
# DT[,Experiment := tolower(Experiment)][str_detect(Proteins,"A0AV96") &Experiment == "g3_gk1_chromatin_eht_3"]


# The HL and ML columns are mixed, where is HL and ML is present then the labels are correct, but if ML is absent from experiment, 
 # the HL is assigned into ML
testing_evidence_sum <- DT[Protein_evidence, on = .(id = id)
                           ]

testing_evidence_sum <- testing_evidence_sum[,.(Prot_HL=median(`Ratio H/L normalized`,na.rm = T)), by = .(majority_protein_i_ds,Experiment)
                             ]
#Why does A0AVT1;A0AVT1-2 from evidence give 0.79601 but in the PG table is 0.799
#are they using all the evidence?
display_venn <- function(x, ...){
    library(VennDiagram)
    grid.newpage()
    venn_object <- venn.diagram(x, filename = NULL, ...)
    grid.draw(venn_object)
}
display_venn( list(
    Protein_evidence_A0AVT1 =  Protein_evidence[majority_protein_i_ds == "A0AVT1;A0AVT1-2",id],
    evidenceDT_evidence_A0AVT1 =  DT[str_detect(Proteins,"A0AVT1"), id]))
#one piece of evidence not in protein evidence?
diff_evidenece <- setdiff(DT[str_detect(Proteins,"A0AVT1"), id],Protein_evidence[majority_protein_i_ds == "A0AVT1;A0AVT1-2",id])
(DT[id == diff_evidenece, ])
#this evidence doesn't belong to this protein group



cols_of_interest <- colnames(Protein_groups_empty) %>% 
    str_subset("Majority|Ratio H/L normalized ")
cols_Type <- colnames(Protein_groups_empty) %>% 
    str_subset("Majority|Ratio H/L normalized |H/L type| H/L iso-count ")


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
Long_protein_groups_type <- Ratio[Type, on = .(majority_protein_i_ds==majority_protein_i_ds, Experiment== Experiment)
][Iso, on = .(majority_protein_i_ds==majority_protein_i_ds, Experiment== Experiment)
][,rm_iso := if_else(iso_count == 0,T,F)]


#  testing_evidence_sum <- testing_evidence_sum[,Experiment := tolower(Experiment)]
# 
# PG_summary_comparison <- testing_evidence_sum[ Long_protein_groups, on = .(Experiment = Experiment,majority_protein_i_ds = majority_protein_i_ds)][
#     , `:=`(Prot_HL = round(Prot_HL,3),
#            Ratio = round(Ratio,3))]
# #to understand disagreement, looking at disagreing proteins with smallest number of evidence
# PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  ][Protein_evidence[,.N, by = .(majority_protein_i_ds)],on = .(majority_protein_i_ds = majority_protein_i_ds)]
# testing_evidence_agreement <- DT[,.(Proteins,`Ratio H/L normalized`, id,Experiment,`Leading Razor Protein`)][Protein_evidence, on = .(id = id)
# ][,Experiment := tolower(Experiment)][majority_protein_i_ds == "P01769" &  Experiment == "g1_px1194_h12" ,]
# #Indeed the median of the two evidences is calculated correctly but still disagrees with the PG summary
# DT[,Experiment := tolower(Experiment)][Proteins == "P01769" &  Experiment == "g1_px1194_h12" ,]
# 
# #This protein has NA in evidence but a value in the PG, why?
# DT[,Experiment := tolower(Experiment)][Proteins == "P01275" &  Experiment == "g1_px309_hcc202_1" ,] #no evidence for this Protein experiment this is because of bad convertion to unique names
# #lets check all the evidence ids of this protein and see if any of them belong to this experiment/ accounted for
# Protein_evidence[majority_protein_i_ds == "P01275",id ]
# DT[id == "2096295" ,] #this shows that the values NaN are not recognised as NA and so never become ML
# 
# 
# testing_evidence_sum <- dcast(testing_evidence_sum, majority_protein_i_ds~Experiment, value.vars = "Prot_HL")
# 
# testing_values <- Protein_groups_empty %>% 
#     janitor::clean_names() %>% 
#     dplyr::select(matches("majority|ratio_h_l_normalized_")) %>% 
#     pivot_longer(cols  = contains("ratio"),names_to = "Experiment", values_to = "Ratio" ) %>% 
#     subset(!is.na(Ratio))
# # Protein_groups_example <- Protein_groups[2,]
# sampling_example <-  sampling_evidence[id %in% str_split(Protein_groups_example$`Evidence IDs`,";",simplify = T) , ]
# 
# sampling_example_protein <-  sampling_evidence[Proteins == Protein_groups_example$`Majority protein IDs` , ]
# sampling_example_protein[,`Ratio M/L normalized`:=as.numeric(`Ratio M/L normalized`)]
# sampling_example_protein_medianed <- sampling_example_protein[, .(median_ratio = median(`Ratio M/L normalized`, na.rm = TRUE)),by = Experiment]

#Dealing with NaNs
# DT <- fread(file_input, select = columns_to_read)
# DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA
# DT <- DT[,`Ratio H/L normalized`:= if_else(is.na(`Ratio H/L normalized`),`Ratio M/L normalized`,`Ratio H/L normalized`)]

testing_evidence_sum <- DT[Protein_evidence, on = .(id = id)
                           ][,leading_majority := str_remove_all(majority_protein_i_ds,";[:graph:]*$")
                             ][`Leading Razor Protein` ==  leading_majority,
                               #][is.finite(`Ratio H/L normalized`)
                                 ]
fixed_exp_names <- data.table(Experiment_old = testing_evidence_sum$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]
testing_evidence_sum <- testing_evidence_sum[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]%>%
    setnames(c("i.Experiment"), c("Experiment"))


testing_evidence_sum_testing <- testing_evidence_sum[majority_protein_i_ds == "A0AV96-2;A0AV96" & Experiment == "g1_gk1_chromatin_al",]
Protein_longer_trial <- Long_protein_groups_type[majority_protein_i_ds == "A0AV96-2;A0AV96" & Experiment == "g1_gk1_chromatin_al",]
testing_evidence_sum_testing[Protein_longer_trial[is.finite(h_l_ratio),], on =.(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)]

testing_evidence_sum <- testing_evidence_sum[Long_protein_groups_type[is.finite(h_l_ratio),], on =.(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)] %>% View()
trial <- testing_evidence_sum[1:100,]
[!(Type == "ISO-MSMS" & rm_iso == T), ]
evidence_n <-testing_evidence_sum[,.N,by = .(Experiment,majority_protein_i_ds)] 
testing_evidence_sum<-testing_evidence_sum[,.(Prot_HL=median(`Ratio H/L normalized`,na.rm = T)), by = .(majority_protein_i_ds,Experiment)
]
evidence_n <- evidence_n[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL] %>%
    setnames(c("i.Experiment"), c("Experiment"))

PG_summary_comparison <- testing_evidence_sum[ Long_protein_groups, on = .(Experiment = Experiment,majority_protein_i_ds = majority_protein_i_ds)][
    , `:=`(Prot_HL = round(Prot_HL,3),
           Ratio = round(Ratio,3))]
#https://groups.google.com/g/maxquant-list/c/06JU00WWfR8/m/SsDGgaWICQAJ
#to understand disagreement, looking at disagreing proteins with smallest number of evidence
PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  
][evidence_n,on = .(majority_protein_i_ds = majority_protein_i_ds,Experiment= Experiment)
][,`:=`(difference = Prot_HL-Ratio)]
plot(log2(PG_summary_comparison$Prot_HL),log2(PG_summary_comparison$Ratio))

PG_summary_comparison[Identical == T,N] %>% .[.<40] %>% hist(breaks = 1000, main = "Hist of evidence number of Correct PG ratio", xlim = c(0,40))
PG_summary_comparison[Identical == F,N] %>% .[.<40] %>% hist(breaks = 1000, main = "Hist of evidence number of Incorrect PG ratio", xlim = c(0,40))

fixed_exp_names <- data.table(Experiment_old = testing_evidence_sum$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]
evidence_n <- evidence_n[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL] %>%
    set_na
testing_evidence_sum <- testing_evidence_sum[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]
testing_evidence_sum <- setnames(testing_evidence_sum, c("i.Experiment"), c("Experiment"))
PG_summary_comparison <- testing_evidence_sum[ Long_protein_groups, on = .(Experiment = Experiment,majority_protein_i_ds = majority_protein_i_ds)][
    , `:=`(Prot_HL = round(Prot_HL,3),
           Ratio = round(Ratio,3))
    ]
Long_protein_groups_type[,h_l_ratio:=NULL]
#https://groups.google.com/g/maxquant-list/c/06JU00WWfR8/m/SsDGgaWICQAJ
#to understand disagreement, looking at disagreing proteins with smallest number of evidence
PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  
                                               ][evidence_n,on = .(majority_protein_i_ds = majority_protein_i_ds,Experiment= Experiment)
                                                 ][,`:=`(difference = Prot_HL-Ratio)
                                                   ][Long_protein_groups_type, on = .(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)
                                                     ][type == "Median",][is.finite(Prot_HL) | is.finite(Ratio),]
experiment = "g3_gk1_chromatin_a_tsa_1"
protein_id = "A1L0T0"
DT[,Experiment := tolower(Experiment)][str_detect(Proteins,protein_id) &  Experiment == experiment ,] 
evidence_for_P01769 <- Protein_evidence[majority_protein_i_ds == protein_id,id]
DT[id %in% evidence_for_P01769 & Experiment == experiment,]#I don't understand what's happening here, the evidence are the same but the calculated values different

DT[,Experiment := tolower(Experiment)][Proteins == "REV__Q9UDX3-2;REV__Q9UDX3;REV__B5MCN3" #&  Experiment == "g1_kw17_130319" 
                                       ,] 
evidence_for_ <- Protein_evidence[majority_protein_i_ds == "REV__Q9UDX3-2;REV__Q9UDX3;REV__B5MCN3",id]
DT[id %in% evidence_for_,]#I don't understand what's happening here, the evidence are the same but the calculated values different


