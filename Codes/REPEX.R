pacman::p_load(tidyverse,data.table,ggpmisc)
#Recreating MQ Protein Intensities
#ProteomeHD exploration#
#Summarise abundance for PGs
#Do peptides that anticorrelate have specific characteristics, could the be put into a "black"list of quantification, do they have more modfications?
#Are all peptide evidence used for quantificaiton?# MULTI-MSMS:     label pair detected
# - MULTI-SECPEP  identification through second peptide algorithm
# - MULTI-MATCH    identification based on match between runs
# - ISO-MSMS:         isotope cluster detected but no label pair
# - MSMS:                no MS1 label state detected
#The type of the feature. 
#'MSMS' - for an MS/MS spectrum without an MS1 isotope pattern assigned.
#ISO-MSMS' - MS1 isotope cluster identified by MS/MS. 
#'MULTI-MSMS' - MS1 labeling cluster identified by MS/MS. 
#MULTI-SECPEP' - MS1 labeling cluster identified by MS/MS as second peptide.
#'MULTI-MATCH' - MS1 labeling cluster identified by matching between runs. 
#'In case of label-free data there is no difference between 'MULTI' and 'ISO'.

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
#it seems that for mot cases, if the HL is NA and the ML is there,
#then the ML is the real HL and when used agrees with the PG ratio
#This is especially true for proteins with a single evidence
DT[,`Ratio H/L normalized`:= fifelse(is.na(`Ratio H/L normalized`) &  (parameter_group != "g3") #& (!is.nan(`Ratio H/L normalized`))
                                     ,`Ratio M/L normalized`,`Ratio H/L normalized`)]

DT[,`Ratio H/L normalized`:= fifelse(is.nan(`Ratio M/L normalized`) &  (parameter_group == "g3") & (Type == "ISO-MSMS"),NaN,`Ratio H/L normalized`)]
Intensity_cols <- DT %>% colnames() %>% str_detect("Intensity ") %>% which()
DT[,Intensity :=fifelse(str_detect(Experiment,"^g3", negate = T), rowSums(.SD, na.rm = T),Intensity),.SDcols =Intensity_cols
        ][,Intensity :=fifelse(Intensity == 0,NaN,Intensity)]
DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA
fixed_exp_names <- data.table(Experiment_old = DT$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]
DT <- DT[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]
DT <- setnames(DT, c("i.Experiment"), c("Experiment"))


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

#
#removing ISOMSMS not used for quantification and counting valid evidence per protein
testing_evidence_sum <- testing_evidence_sum[Long_protein_groups_type[is.finite(h_l_ratio),], on = .(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)]
finding_ISOMS_rule <- testing_evidence_sum
testing_evidence_sum <- testing_evidence_sum[!(rm_iso == T & Type == "ISO-MSMS")][!is.na(`Ratio H/L normalized`)]
evidence_n <-testing_evidence_sum[,.N,by = .(Experiment,majority_protein_i_ds)] 

#median of valid evidence
testing_evidence_sum<-testing_evidence_sum[,.(Prot_HL=median(log2(`Ratio H/L normalized`),na.rm = T)), by = .(majority_protein_i_ds,Experiment)
]

#combine with PG ratio
PG_summary_comparison <- testing_evidence_sum[ Long_protein_groups, on = .(Experiment = Experiment,majority_protein_i_ds = majority_protein_i_ds)][
  , `:=`(Prot_HL = round(Prot_HL,2),
         Ratio = round(log2(Ratio),2))
]

#removing proteins which didn't arise from median (advanced ratio calculation)
PG_summary_comparison <- PG_summary_comparison[,Identical := .(Prot_HL==Ratio)  
][evidence_n,on = .(majority_protein_i_ds = majority_protein_i_ds,Experiment= Experiment)
][,`:=`(difference = Prot_HL-Ratio)
][Long_protein_groups_type, on = .(majority_protein_i_ds = majority_protein_i_ds, Experiment = Experiment)
][!(is.na(Prot_HL) &is.na(Ratio)),
  ][,parameter_group:=str_match(Experiment,"^([:graph:]*?)_") %>% .[,2]]  
testing <- PG_summary_comparison[N<10  & parameter_group =="g3",.(N, Prot_HL, Ratio,type)] %>% 
  mutate(N = as.character(N)  ) 
my.formula <- y ~ x
p <- ggplot(data = testing, aes(x = Prot_HL, y = Ratio, colour = type)) +
  # geom_smooth(method = "lm", se=FALSE, formula = my.formula, alpha = 0.2) +
  stat_poly_eq(formula = my.formula,size= 3,
               aes(label = paste(..rr.label..)), 
               parse = TRUE) +         
  geom_point(alpha = 0.2)+
  facet_wrap("N")+
  lims(x= c(-13,14), y = c(-15,15))
p

#plots
PG_summary_comparison[type == "Median",][Identical == T,N] %>% .[.<40] %>% 
  hist(breaks = 1000, main = "Hist of evidence number of Correct PG ratio", xlim = c(0,40))
PG_summary_comparison[type == "Median",][Identical == F,N] %>% .[.<40] %>% hist(breaks = 1000, main = "Hist of evidence number of Incorrect PG ratio", xlim = c(0,40))

PG_summary_comparison[!(is.na(Prot_HL)|is.na(Ratio))] %>%  
  ggplot(aes(x = Prot_HL, y = Ratio, colour= fifelse(N%%2 == 0,"Even","Odd")))+
  geom_point(alpha = 0.3)+labs(title="Recalculating ProteinGroup Silac Ratio from Evidence file",
                               x ="Ratio Recalculated from Evidence file", y = "Ratio from ProteinGroup",colour = "Number of Evidence for PG in Experiment")

ggsave(here::here("Output","Silac_Ratio_PG_reconstruction.png"))
PG_summary_comparison[,parameter_group:= Experiment %>% str_match("^([:graph:]*?)_") %>% .[,2]]

PG_summary_comparison[N==2,][!(is.na(Prot_HL)|is.na(Ratio))] %>%  
  ggplot(aes(x = Prot_HL, y = Ratio, colour= parameter_group))+
  
  geom_point(alpha = 0.3)+facet_wrap(~parameter_group)+labs(title="N == 2 based on parameter group",
                                                            x ="Ratio Recalculated from Evidence file", y = "Ratio from ProteinGroup",colour = "parameter group")

DT[,Experiment := tolower(Experiment)]
DT_imputted <- DT   
experiment = "g3_kw35_ne_2"
protein_id = "^O14950"


DT[str_detect(Proteins,protein_id) &  str_detect(Experiment,experiment)# & Type != "ISO-MSMS"
     ,`Ratio H/L normalized`] 
`DT_imputted[str_detect(Proteins,protein_id) &  str_detect(Experiment,experiment) & Type != "ISO-MSMS"
   ,] 
 testing_evidence_sum[str_detect(majority_protein_i_ds,protein_id) &  str_detect(Experiment,experiment)]
 evidence_ids <- DT[str_detect(Proteins,protein_id) &  str_detect(Experiment,experiment),id]
 Protein_evidence[id %chin% evidence_ids]

PG_summary_comparison[N==2,]
####Advanced  Ration Estimation
# columns_to_read <- c("Sequence", "Proteins", "Experiment",# "Protein group IDs",
#                      "Leading Razor Protein","id","Intensity",
#                      "Ratio H/L normalized", "Ratio M/L normalized",# "Ratio H/M normalized" ,"Peptide ID",
#                      # "Reverse", "Potential contaminant",
#                      "Type")

to_add_DT <- DT[is.finite(Intensity) & is.finite(`Ratio H/L normalized`)]#[parameter_group == "g3"]
Regression_df <-to_add_DT[Long_protein_groups_type[is.finite(h_l_ratio),.(majority_protein_i_ds,Experiment, rm_iso)
                                             ], 
                    on = .(Proteins = majority_protein_i_ds, Experiment = Experiment)
                    ][!(Type == "ISO-MSMS" & rm_iso == T)][, n:=.N, by=.(Proteins,Experiment), 
                                       ][is.finite(Intensity) & is.finite(`Ratio H/L normalized`)][ n > 3]
Regression_ <- Regression_df[, {LM = lm(log2(Intensity)~log2(`Ratio H/L normalized`), data = .SD)
LM.summary = summary(LM)
list(Intercept = LM$coefficients[1],
     RatioHL = LM$coefficients[2],
     R2 = LM.summary$r.squared,
     pvalue.Inter = LM.summary$coefficients["(Intercept)", 4],
     pvalue.RatioHL = LM.summary$coefficients["log2(`Ratio H/L normalized`)", 4])
},
by = .(Proteins,Experiment)]
Regression_types <- Regression_[Long_protein_groups_type[is.finite(h_l_ratio),.(majority_protein_i_ds,Experiment, type)
],on = .(Proteins = majority_protein_i_ds, Experiment = Experiment)][is.finite(Intercept)]


  ggplot(Regression_types[,parameter_group:=Experiment %>% str_match("^([:graph:]*?)_") %>% .[,2]][
    between(RatioHL,-15,15 )],#[,head(.SD,2541),by = .(type,parameter_group)]
         aes(x= R2, 
             y = -log10(pvalue.RatioHL),
             colour = RatioHL))+
    geom_point(alpha = 0.5)+
    scale_color_gradient2(mid = "grey95")+
    theme_bw()+
    facet_wrap(type~parameter_group, scales = "free_y")
N2_DT <- melt(Regres sion_types[,#`pvalue.RatioHL`<0.05 & RatioHL>0
                                colnames(Regression_types), with=FALSE], id.vars = c("Proteins","Experiment","type"))
N2_DT[,value := fifelse(variable == "pvalue.Inter", log10(value),value)]
ggplot(N2_DT[is.finite(value) &
               !(variable == "RatioHL" & !between(value,-15,15 ))&
               !(variable == "Intercept" & !between(value,-2,30 ))],aes(x = value))+
  geom_density(aes(fill=type),alpha = 0.5, position = 'identity', bins = 15)+
  facet_wrap(~variable,scales = "free")+ 
  # scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    ggtitle("Linear Regression output as features for Decision tree for Adv. Rat. Est")

#Setup a binary classification problem
pacman::p_load(DAAG,party,rpart,rpart.plot,mlbench,caret,pROC,tree)
mydata <- Regression_types[type != "Plateau",.(Intercept,RatioHL,R2,pvalue.Inter,pvalue.RatioHL,type)
                           ][,head(.SD,2541),by = type]
set.seed(1234)
ind <- sample(2, nrow(mydata), replace = T, prob = c(0.75,0.25))
train <- mydata[ind == 1,]
test <- mydata[ind == 2,]
tree <- rpart(type ~., data = train)
rpart.plot(tree)

p1 <- predict(tree, test, type = 'prob')
p1 <- p1[,2]
r <- multiclass.roc(test$type, p1, percent = TRUE)
roc <- r[['rocs']]
r1 <- roc[[1]]
plot.roc(r1,
         print.auc=TRUE,
         auc.polygon=TRUE,
         grid=c(0.1, 0.2),
         grid.col=c("green", "red"),
         max.auc.polygon=TRUE,
         auc.polygon.col="lightblue",
         print.thres=TRUE,
         main= 'ROC Curve')

p <- predict(tree, train, type = 'class')
confusionMatrix(p, train$type %>% as.factor(), positive = "Median")

#Figuring ISOMSMS RM 
library(e1071)
ISO_MS_N_Ratio <-  finding_ISOMS_rule[,.(`Leading Razor Protein`,
                      Type,Experiment,`Ratio H/L normalized`,iso_count)
                   ][,is_ISOMSMS := fifelse(Type == "ISO-MSMS",T,F)
                     ][,sum_ISO := sum(is_ISOMSMS), by = .(`Leading Razor Protein`, Experiment)
                       ][sum_ISO>0][,evidence_n := .N, by = .(`Leading Razor Protein`, Experiment)
                       ][,.(`Leading Razor Protein`,evidence_n,sum_ISO,
                          Experiment,iso_count)] %>% unique() %>% 
  .[,Ratio_iso_total := sum_ISO/evidence_n]

ISO_MS_N_Ratio[,iso_used := fifelse(iso_count>0,T,F)]

dat = ISO_MS_N_Ratio[,.(Ratio_iso_total,evidence_n,iso_used)
                     ][,Ratio_iso_total :=Ratio_iso_total^3
                       ][,iso_used:= as.factor(iso_used)][,head(.SD,230004),by = iso_used]
# Splitting the dataset into the Training set and Test set
#install.packages('caTools')
library(caTools)

set.seed(123)
split = sample.split(dat$iso_used, SplitRatio = 0.75)

training_set = subset(dat, split == TRUE)
test_set = subset(dat, split == FALSE)

classifier = e1071::svm(formula = iso_used ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'radial')
y_pred = predict(classifier, newdata = test_set[,.(evidence_n,Ratio_iso_total)])
save(classifier,file = here::here("Datasets","Processed","ISO_SVM_ProteomeHD.rda"))
cm = table(test_set$iso_used , y_pred)
cm

roc_svm_test <- pROC::roc(response = as.numeric(test_set$iso_used), predictor =as.numeric(y_pred))
plot.new()
plot(roc_svm_test, add = TRUE,col = "red", print.auc=TRUE, print.auc.x = 0.5, print.auc.y = 0.3)
legend(0.3, 0.2, legend = c("test-svm"), lty = c(1), col = c("blue"))

ISO_MS_N_Ratio %>% 
  ggplot(aes(x = Ratio_iso_total, colour = iso_used, y = log2(evidence_n)))+
  geom_point(alpha = 0.5)+
  ggtitle("Amount of evidence needed to discard ISO-MSMS")


