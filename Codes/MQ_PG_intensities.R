
columns_to_read <- c("Sequence", "Proteins", "Experiment",# "Protein group IDs",
                     "Leading Razor Protein","id",
                     "Ratio H/L normalized", "Ratio M/L normalized",# "Ratio H/M normalized" ,"Peptide ID",
                     # "Reverse", "Potential contaminant",
                     "Type")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")

# Read in data
DT <- fread(file_input, select = columns_to_read)
DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA
DT[,`Ratio H/L normalized`:= if_else(is.na(`Ratio H/L normalized`) & !is.nan(`Ratio H/L normalized`),`Ratio M/L normalized`,`Ratio H/L normalized`)]
DT[,`Ratio H/L normalized`:= if_else(is.nan(`Ratio M/L normalized`),NaN,`Ratio H/L normalized`)]


DT_cols <- fread(file_input, nrows = 1)

Protein_groups_empty <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), nrows = 2)
Protein_groups <- fread(input = here::here("Datasets","Raw","ProteomeHD", "proteinGroups.txt"), select = c("Intensity", "Protein IDs","Majority protein IDs", "Evidence IDs", "Peptide IDs"))

Protein_evidence <- Protein_groups %>%
    janitor::clean_names() %>% 
    dplyr::select(c(majority_protein_i_ds,evidence_i_ds)) %>% 
    separate_rows(evidence_i_ds) %>% 
    rename(id = evidence_i_ds) %>% 
    setDT() 
Protein_evidence[,id := as.numeric(id)]

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




testing_evidence_sum <- DT[Protein_evidence, on = .(id = id)
][,leading_majority := str_remove_all(majority_protein_i_ds,";[:graph:]*$")
][`Leading Razor Protein` ==  leading_majority,
  ][is.finite(`Ratio H/L normalized`)
]
fixed_exp_names <- data.table(Experiment_old = testing_evidence_sum$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]

testing_evidence_sum <- testing_evidence_sum[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]%>%
    setnames(c("i.Experiment"), c("Experiment"))

testing_evidence_sum <- testing_evidence_sum[Long_protein_groups_type[is.finite(h_l_ratio),], 
                                             on =.(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)]
testing_evidence_sum <- testing_evidence_sum[!(Type == "ISO-MSMS" & rm_iso == T), ][type == "Median",]
evidence_n <-testing_evidence_sum[,.N,by = .(Experiment,majority_protein_i_ds)] 
testing_evidence_sum <- testing_evidence_sum[,.(Prot_HL=median(`Ratio H/L normalized`,na.rm = T)), by = .(majority_protein_i_ds,Experiment)
]

PG_summary_comparison <- testing_evidence_sum[ Long_protein_groups, on = .(Experiment = Experiment,majority_protein_i_ds = majority_protein_i_ds)][
    , `:=`(Prot_HL = round(Prot_HL,3),
           Ratio = round(Ratio,3))]
#https://groups.google.com/g/maxquant-list/c/06JU00WWfR8/m/SsDGgaWICQAJ
#to understand disagreement, looking at disagreing proteins with smallest number of evidence
PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  
][evidence_n,on = .(majority_protein_i_ds = majority_protein_i_ds,Experiment= Experiment)
][,`:=`(difference = Prot_HL-Ratio,
        Even_evidence = if_else(N%%2 ==1,"odd","even"))]
ggplot(PG_summary_comparison,
       aes(x = log2(Prot_HL), y = log2(Ratio), colour = Even_evidence))+
  geom_point()+
  ggtitle("ProteomeHD recalculated PG ratios",
          subtitle = "Linear and Plateau removed, even evidence still an issue")
plot(log2(PG_summary_comparison$Prot_HL),log2(PG_summary_comparison$Ratio))

PG_summary_comparison[Identical == T,] %>% pull(N) %>% .[.<40] %>% hist(breaks = 1000, main = "Hist of evidence number of Correct PG ratio", xlim = c(0,40))
PG_summary_comparison[Identical == F,] %>% pull(N)%>% .[.<40] %>% hist(breaks = 1000, main = "Hist of evidence number of Incorrect PG ratio", xlim = c(0,40))

Long_protein_groups_type[,h_l_ratio:=NULL]
#https://groups.google.com/g/maxquant-list/c/06JU00WWfR8/m/SsDGgaWICQAJ
#to understand disagreement, looking at disagreing proteins with smallest number of evidence
PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  
][evidence_n,on = .(majority_protein_i_ds = majority_protein_i_ds,Experiment= Experiment)
][,`:=`(difference = Prot_HL-Ratio)
][Long_protein_groups_type, on = .(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)
][type == "Median",][is.finite(Prot_HL) | is.finite(Ratio),]


