####Protein_Peptide_Cor_Heatmaps####
pacman::p_load(tidyverse,data.table,seqinr)
Protein_Sequences <- read.fasta(file = here::here("Datasets","Raw","ProteomeHD", "Uniprot_human_GK_ProtHD_280515.fasta"), 
                                seqtype = "AA", as.string = TRUE) %>% 
    unlist %>% enframe(name = "Uniprot",value = "Prot_Sequence") %>% as.data.table %>% 
    setkey("Uniprot")
Protein_Sequences[,Uniprot:= str_match(Uniprot,"\\|([:graph:]*?)\\|")[,2]]
columns_to_read <- c("Sequence", "Proteins", "Experiment",
                     "Leading Razor Protein","id","Intensity","Intensity L","Intensity M", "Intensity H",
                     "Ratio H/L normalized", "Ratio M/L normalized","Type")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")
columns_peptides <- c("Sequence", "Proteins",
                      "Leading razor protein","Start position","End position")

# file_peptides <- here::here("Datasets","Raw","ProteomeHD", "peptides.txt")
# Peptides_cols <- fread(file_peptides, select = columns_peptides)

# Read in data
DT <- fread(file_input, select = columns_to_read)

#The NA meaning changes based on double or triple silac
DT[,parameter_group:= Experiment %>% str_match("^([:graph:]*?)_") %>% .[,2]]

#ML ratio in double silac is actualy HL
DT[,`Ratio H/L normalized`:= fifelse(is.na(`Ratio H/L normalized`) &  (parameter_group != "g3") #& (!is.nan(`Ratio H/L normalized`))
                                     ,`Ratio M/L normalized`,`Ratio H/L normalized`)]
DT[,`Ratio H/L normalized`:= fifelse(is.nan(`Ratio M/L normalized`) &  (parameter_group == "g3") & 
                                         (Type == "ISO-MSMS"),NaN,`Ratio H/L normalized`)]
##

#This needs to also be fixed on the Intensity Cols
Intensity_cols <- DT %>% colnames() %>% str_detect("Intensity ") %>% which()
DT[,Intensity :=fifelse(str_detect(Experiment,"^g3", negate = T), rowSums(.SD, na.rm = T),Intensity),.SDcols =Intensity_cols
][,Intensity :=fifelse(Intensity == 0,NaN,Intensity)]
##

DT$`Ratio H/L normalized`[!is.finite(DT$`Ratio H/L normalized`)] <- NA #Convert NaN to NA

#make experiment names consistent with lowercase and no gaps
fixed_exp_names <- data.table(Experiment_old = DT$Experiment %>% unique)[,Experiment := janitor::make_clean_names(Experiment_old)]
DT <- DT[fixed_exp_names, on = .(Experiment= Experiment_old)][,Experiment:=NULL]
DT <- setnames(DT, c("i.Experiment"), c("Experiment"))
#

# For each peptide, get median M/L ratio per Experiment
med_DT <- DT[, .(median_ratio = median(`Ratio H/L normalized`, na.rm = TRUE)),
                            .(Sequence, Proteins, Experiment, `Leading Razor Protein`)]
med_DT[,n:=.N, by = c("Proteins","Experiment")]
#med_DT <- med_DT[n>10] #too few evidence are removed

# Exlude peptides that were seen in less than 25 experiments, since we can't create a solid correlation profile
exclude_peptides <- med_DT[, .N, Sequence][ N <= 25, Sequence ]
med_DT <- med_DT[ !Sequence %in% exclude_peptides][, median_ratio := log2(median_ratio)]
med_DT[,distinct_groups := uniqueN(Proteins), by = .(`Leading Razor Protein`)] #how many PG in the same Leading PG
# med_DT <- med_DT[distinct_groups>1]

# Calculate correlations among peptides for each protein
pb <- txtProgressBar(min = 0, max =length(unique(med_DT$`Leading Razor Protein`)), style = 3)
correlations <- data.table()
#we should be cycling through Uniprot Ids, not Protein Groups
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))
Proteins_of_interest <- "O00303,O15371,O15372,O75821,O75822,P55884-2,P60228,Q13347,Q14152,Q5JSZ5,Q7L2H7,Q99613,Q9UBQ5,Q9Y262" %>% 
    str_split(",") %>% unlist()
save_pheatmap_pdf <- function(x, filename, width=480, height=480) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    png(filename, width=width, height=height,units = "px")
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

for(i in  Proteins_of_interest){
    tmp_protein <-  i# unique(med_DT$`Leading Razor Protein`)[i] #""#  #
    tmp <- med_DT[str_detect(Proteins,tmp_protein)]#med_DT[`Leading Razor Protein` == tmp_protein]
 
    tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
    tmp$Experiment <- NULL
    #finding and removing correlations with less than 10 common datapoints
    missingness_tmp <- crossprod(is.na(tmp))>10
    tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")#using suggested bicor
    tmp[missingness_tmp == FALSE] <- 0
    #tmp[is.na(tmp)] <- 0
    is_0 <- is.na(tmp)
    keep_few_0 <- is_0 %>% matrixStats::rowSums2() %>% set_names(rownames(tmp)) %>% .[.<(ncol(tmp)/1.5)] %>% names() #only keep peptides with few 0/NA
    if(length(keep_few_0)>2){
        tmp_heatmap <- tmp[keep_few_0,keep_few_0]
        tmp_prot_seq <- Protein_Sequences[Uniprot == tmp_protein,Prot_Sequence]
        tmp_pep_col <- DT[str_detect(Proteins,tmp_protein),.(Sequence,Proteins,`Leading Razor Protein`) 
        ]%>% unique() %>% .[,Midpoint:= map_dbl(.x = Sequence,~str_locate(tmp_prot_seq,.x) %>% mean())] 
        
        
        
        annotation <- tmp_pep_col[,.(Sequence,Proteins,`Leading Razor Protein`, Midpoint)] %>% column_to_rownames("Sequence") %>% 
            set_names(c("PG","Majority","Peptide Position")) %>% 
            mutate(`Peptide Position` = rank(`Peptide Position`))
        protein_order <- annotation %>% subset(rownames(.) %in% rownames(tmp_heatmap)) %>% arrange(`Peptide Position`) %>% rownames()
        
       xx <-  pheatmap::pheatmap(tmp_heatmap[protein_order,protein_order],show_rownames = F,show_colnames = F,
                           cluster_rows = F,cluster_cols = F,
                           annotation = annotation, color=myColor, breaks=myBreaks, main = tmp_protein)
       save_pheatmap_pdf(xx, glue::glue("{tmp_protein}_pep_cor.png"))
    }
    
    tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))
    tmp <- tmp[, .( Peptide1 = as.character(Var1), Peptide2 = as.character(Var2), PCC = value ) ]
    tmp <- tmp[ Peptide1 >= Peptide2 ][!is.nan(PCC)] #%>% distinct()
    tmp[, Protein := tmp_protein ]
    correlations <- rbind( correlations, tmp)
    setTxtProgressBar(pb, i)
}
