#MaxQuant Peptide Inference##
#Looking at the evidence file
#Peptides can be assigned to a single or multiple proteins
#And proteins in proteins groups can individually exist for other peptides
# sp|P62874|GBB1_MOUSE;sp|P62880|GBB2_MOUSE is an interesting case where MaxQuant assigns them to GBB2 even though they
# have the same number of entries, with GBB1 having more unique
# rawfile level and project level inference is different
# quantitation3
evidence_original <- janitor::clean_names(evidence)
evidence <- evidence_original %>% 
    subset(!duplicated(sequence))
all_shared <- evidence %>% 
    pull(proteins) %>%     unique() %>% 
    str_subset(";") 
example <-  all_shared[5]
maxquant_razor <- evidence %>% 
    subset(str_detect(proteins, ";") & !duplicated(proteins)) 
    

calculated_Protein_peptides <- function(protein_Group_entry){
    # protein_Group_entry <-   "sp|P62715|PP2AB_MOUSE;sp|P63330|PP2AA_MOUSE"
    split_protein_Group_entry <- str_split(protein_Group_entry,";", simplify = T) %>% 
        str_replace_all("\\|","\\\\|") %>% set_names(str_split(protein_Group_entry,";", simplify = T))
    evidence_for_protein <- purrr::map_dbl(.x = split_protein_Group_entry,
               ~subset(evidence,str_detect(proteins,.x) & (proteins != protein_Group_entry)) %>% nrow()) %>% 
        sort(decreasing = T)
    # data.frame(proteins = protein_Group_entry,
    #            leading_razor_protein = evidence_for_protein)
    data.frame(proteins = protein_Group_entry,
               leading_razor_protein = names(evidence_for_protein[1]))
}
testing <- purrr::map_dfr(maxquant_razor$proteins %>% 
                              unique(),calculated_Protein_peptides)
testing_confict <- purrr::map_dfr(subset(maxquant_razor, 
                                         maxquant_razor$leading_razor_protein != testing$leading_razor_protein) %>% 
                                      pull(proteins),calculated_Protein_peptides)

(maxquant_razor$leading_razor_protein==testing$leading_razor_protein) %>% table()
