##### evıdence.txt to proteın groups.txt ####
load(here::here("Datasets","Processed","ISO_SVM_ProteomeHD.rda"))
columns_to_read <- c("Sequence", "Proteins", "Experiment",# "Protein group IDs",
                     "Leading Razor Protein","id","Intensity","Intensity L","Intensity M", "Intensity H",
                     "Ratio H/L normalized", "Ratio M/L normalized",# "Ratio H/M normalized" ,"Peptide ID",
                     # "Reverse", "Potential contaminant",
                     "Type")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")

# Read in data
DT <- fread(file_input, select = columns_to_read)
#DT <- DT[str_detect(Experiment,"^g3_"),]
DT[,parameter_group:= Experiment %>% str_match("^([:graph:]*?)_") %>% .[,2]]
if(str_detect(file_input,"ProteomeHD")){
DT[,`Ratio H/L normalized`:= fifelse(is.na(`Ratio H/L normalized`) &  (parameter_group != "g3") #& (!is.nan(`Ratio H/L normalized`))
                                     ,`Ratio M/L normalized`,`Ratio H/L normalized`)]

DT[,`Ratio H/L normalized`:= fifelse(is.nan(`Ratio M/L normalized`) &  (parameter_group == "g3") & (Type == "ISO-MSMS"),NaN,`Ratio H/L normalized`)]
Intensity_cols <- DT %>% colnames() %>% str_detect("Intensity ") %>% which()
DT[,Intensity :=fifelse(str_detect(Experiment,"^g3", negate = T), rowSums(.SD, na.rm = T),Intensity),.SDcols =Intensity_cols
][,Intensity :=fifelse(Intensity == 0,NaN,Intensity)]
DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA
}
fixed_exp_names <- data.table(Experiment_old = DT$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]
DT <- DT[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]
DT <- setnames(DT, c("i.Experiment"), c("Experiment"))


#removing entries which are not needed because of ISOMSMS and enough other evidence
Predict_ISOMS <- DT[,.(Type,Experiment,`Leading Razor Protein`)][,is_ISOMSMS := fifelse(Type == "ISO-MSMS",T,F)
                                               ][,sum_ISO := sum(is_ISOMSMS), by = .(`Leading Razor Protein`, Experiment)
                                                 ][,evidence_n := .N, by = .(`Leading Razor Protein`, Experiment)
                                                   ] %>% unique() %>% 
       .[,Ratio_iso_total := sum_ISO/evidence_n] 
Predict_ISOMS$Use_ISO = predict(classifier, newdata = Predict_ISOMS[,.(evidence_n,Ratio_iso_total)])
Predict_ISOMS <- Predict_ISOMS[,.(Use_ISO,`Leading Razor Protein`,Experiment)] %>% unique()
DT <- Predict_ISOMS[DT, on = .(`Leading Razor Protein` = `Leading Razor Protein`, Experiment = Experiment)]
DT <- DT[!(Type == "ISO-MSMS" & Use_ISO == F)]

PG_Ratio_HL <- DT[is.finite(`Ratio H/L normalized`)
                           ][,.(Prot_HL=median(log2(`Ratio H/L normalized`),na.rm = T)), by = .(`Leading Razor Protein`,Experiment)
                             ]
