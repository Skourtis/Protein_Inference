####Correlating peptide evidence####
#Should we use sequence level median, before PG median? - agrees with Peptide correlation median
#How do we select peptides to classify as a more or less trustworthy
#what do we do with odd and event
pacman::p_load(tidyverse,data.table,
               # ,amap
               seqinr)
Protein_Sequences <- read.fasta(file = here::here("Datasets","Raw","human.fasta"), 
           seqtype = "AA", as.string = TRUE, forceDNAtolower = TRUE,
           set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE,
           whole.header = FALSE,
           bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
           endian = .Platform$endian, apply.mask = TRUE) %>% 
    unlist %>% enframe(name = "Uniprot",value = "Prot_Sequence") %>% as.data.table %>% 
    setkey("Uniprot")
Protein_Sequences[,Uniprot:= str_match(Uniprot,"\\|([:graph:]*?)\\|")[,2]]
columns_to_read <- c("Sequence", "Proteins", "Experiment",# "Protein group IDs",
                     "Leading Razor Protein","id","Intensity","Intensity L","Intensity M", "Intensity H",
                     "Ratio H/L normalized", "Ratio M/L normalized",# "Ratio H/M normalized" ,"Peptide ID",
                     # "Reverse", "Potential contaminant",
                     "Type")
file_input <- here::here("Datasets","Raw","ProteomeHD", "evidence.txt")
columns_peptides <- c("Sequence", "Proteins",
                     "Leading razor protein","Start position","End position")

file_peptides <- here::here("Datasets","Raw","ProteomeHD", "peptides.txt")
Peptides_cols <- fread(file_peptides, select = columns_peptides)[,Midpoint:= median(c(`Start position`,`End position`)),by = Sequence]

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

multiple_isoforms <- DT[str_detect(Proteins,";") ] #7021 leading protein groups with multiple isoforms

# For each peptide, get median M/L ratio per Experiment
med_DT <- multiple_isoforms[, .(median_ratio = median(`Ratio H/L normalized`, na.rm = TRUE)),
                            .(Sequence, Proteins, Experiment, `Leading Razor Protein`)]
med_DT[,n:=.N, by = c("Proteins","Experiment")]
# med_DT <- med_DT[n>10] #733 leading razor with more than 20
pep_PG <- unique(med_DT[,.(Sequence, Proteins, `Leading Razor Protein`)])
pep_PG <- pep_PG[,unique:= if_else( Proteins ==`Leading Razor Protein`,"unique","non_unique" )][unique == "unique",.(Sequence)][["Sequence"]]

# Exlude peptides that were seen in less than 25 experiments, since we can't create a solid correlation profile
exclude_peptides <- med_DT[, .N, Sequence][ N <= 25, Sequence ]
med_DT <- med_DT[ !Sequence %in% exclude_peptides][, median_ratio := log2(median_ratio)]
med_DT[,distinct_groups := uniqueN(Proteins), by = .(`Leading Razor Protein`)]
# med_DT <- med_DT[distinct_groups>1]
# Calculate correlations among peptides for each protein
pb <- txtProgressBar(min = 0, max =length(unique(med_DT$`Leading Razor Protein`)), style = 3)
correlations <- data.table()
#maybe instead of cycling through Proteins, we should cycle through Leading Razor proteins,
#Because the Proteins and LEading proteins are treated differently here but will be used for the quant of the same Protein
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))

for(i in  1:length(unique(med_DT$`Leading Razor Protein`))){
     i <- 23
    tmp_protein <-  "Q9UPN3"# unique(med_DT$`Leading Razor Protein`)[i] #""#  #
     tmp <- med_DT[`Leading Razor Protein` == tmp_protein]
    # tmp <- med_DT[ `Leading Razor Protein` == tmp_protein ]
     # dcast(tmp,Experiment ~ Sequence, value.var = "median_ratio" ) %>% melt(id.vars = "Experiment", measure.vars = colnames(.)[-1]) %>% 
     #     .[annotation,on = .(variable= Sequence)] %>% 
     #     ggplot(aes(x  = Experiment, y = value, group = variable, colour = Cluster))+geom_point()+
     #     geom_line()+
     #     ggtitle(tmp_protein)
    tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
    tmp$Experiment <- NULL
     #finding and removing correlations with less than 10 common datapoints
    missingness_tmp <- crossprod(is.na(tmp))>10
    # Distance <- amap::Dist(t(tmp),method = "correlation",diag = TRUE, upper = TRUE)
    # Distance[missingness_tmp == FALSE] <- 1
    # Distance[is.na(Distance)] <- 1
    # # Distance[is.na(Distance)] <- 1
    tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")#using suggested bicor
    tmp[missingness_tmp == FALSE] <- 0
    # tmp[is.na(tmp)] <- 0
    is_0 <- is.na(tmp)
    keep_few_0 <- is_0 %>% matrixStats::rowSums2() %>% set_names(rownames(tmp)) %>% .[.<(ncol(tmp)/1.5)] %>% names() #only keep peptides with few 0/NA
    if(length(keep_few_0)>2){
    tmp_heatmap <- tmp[keep_few_0,keep_few_0]
    #   res <- dbscan::optics(as.dist(1-tmp), minPts = 5)
    #  res <- dbscan::extractXi(res, xi = 0.05)
    # res$clusters_xi
    # 
    # hdbscan_peptides <- dbscan::hdbscan(Distance, minPts = 5)
    # annotation <- data.frame(Cluster = res$cluster
    #                              # hdbscan_peptides$cluster
    #                          ,
    #                          Seq = colnames(tmp) )%>% column_to_rownames("Seq")
    # names_order <- annotation %>% arrange(Cluster) %>% rownames()
    tmp_prot_seq <- Protein_Sequences[Uniprot == "Q9UPN3",Prot_Sequence]
    tmp_pep_col <- Peptides_cols[str_detect(Proteins,tmp_protein),.(Sequence,Proteins,`Leading razor protein`)
                                 ][,Midpoint:= map_dbl(.x = Sequence,~str_locate(tmp_prot_seq,.x) %>% mean())] 
    
    

    annotation <- tmp_pep_col[,.(Sequence,Proteins,`Leading razor protein`, Midpoint)] %>% column_to_rownames("Sequence") %>% 
        set_names(c("PG","Majority","Peptide Position")) %>% 
        mutate(
            # Isoforms = case_when(
            # str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4")~"Q9UPN3-5;Q9UPN3-4",
            # str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4",negate = T)~"Q9UPN3-5",
            # str_detect(Isoforms,"-5",negate = T) & str_detect(Isoforms,"-4")~"Q9UPN3-4",        ),
               `Peptide Position` = rank(`Peptide Position`))
    protein_order <- annotation %>% subset(rownames(.) %in% rownames(tmp_heatmap)) %>% arrange(`Peptide Position`) %>% rownames()
    # tmp <- tmp %>% as.matrix %>% .[names_order,names_order]
    # tmp[is.na(tmp)] <- 0
    pheatmap::pheatmap(tmp_heatmap[protein_order,protein_order],show_rownames = F,show_colnames = F,
                        cluster_rows = F,cluster_cols = F,
                       annotation = annotation, color=myColor, breaks=myBreaks, main = tmp_protein)
    }
    
    # testing <- amap::Dist(tmp,"correlation")
    
    

    tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))
    tmp <- tmp[, .( Peptide1 = as.character(Var1), Peptide2 = as.character(Var2), PCC = value ) ]
    tmp <- tmp[ Peptide1 >= Peptide2 ][!is.nan(PCC)] #%>% distinct()
    #if using leading razor, then we can then look back at a simplified evidence table to see if the peptides were unique or shared in Proteins
    tmp[, Protein := tmp_protein ]
    correlations <- rbind( correlations, tmp)
    setTxtProgressBar(pb, i)
}
for(i in   4000:length(unique(med_DT$`Leading Razor Protein`))){
    # i <- 23
    tmp_protein <-  "P01860" #unique(med_DT$`Leading Razor Protein`)[i] #""#  #"Q9UPN3"# #
    tmp <- med_DT[str_detect(Proteins,tmp_protein)]
    # tmp <- med_DT[ `Leading Razor Protein` == tmp_protein ]
    # dcast(tmp,Experiment ~ Sequence, value.var = "median_ratio" ) %>% melt(id.vars = "Experiment", measure.vars = colnames(.)[-1]) %>% 
    #     .[annotation,on = .(variable= Sequence)] %>% 
    #     ggplot(aes(x  = Experiment, y = value, group = variable, colour = Cluster))+geom_point()+
    #     geom_line()+
    #     ggtitle(tmp_protein)
    tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
    tmp$Experiment <- NULL
    #finding and removing correlations with less than 10 common datapoints
    missingness_tmp <- crossprod(is.na(tmp))>10
    # Distance <- amap::Dist(t(tmp),method = "correlation",diag = TRUE, upper = TRUE)
    # Distance[missingness_tmp == FALSE] <- 1
    # Distance[is.na(Distance)] <- 1
    # # Distance[is.na(Distance)] <- 1
    tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")#using suggested bicor
    tmp[missingness_tmp == FALSE] <- 0
    # tmp[is.na(tmp)] <- 0
    is_0 <- is.na(tmp)
    keep_few_0 <- is_0 %>% matrixStats::rowSums2() %>% set_names(rownames(tmp)) %>% .[.<(ncol(tmp)/1.5)] %>% names() #only keep peptides with few 0/NA
    if(length(keep_few_0)>2){
        tmp_heatmap <- tmp[keep_few_0,keep_few_0]
        #   res <- dbscan::optics(as.dist(1-tmp), minPts = 5)
        #  res <- dbscan::extractXi(res, xi = 0.05)
        # res$clusters_xi
        # 
        # hdbscan_peptides <- dbscan::hdbscan(Distance, minPts = 5)
        # annotation <- data.frame(Cluster = res$cluster
        #                              # hdbscan_peptides$cluster
        #                          ,
        #                          Seq = colnames(tmp) )%>% column_to_rownames("Seq")
        # names_order <- annotation %>% arrange(Cluster) %>% rownames()
        tmp_prot_seq <- Protein_Sequences[Uniprot == tmp_protein,Prot_Sequence]
        tmp_pep_col <- Peptides_cols[str_detect(Proteins,tmp_protein),.(Sequence,Proteins,`Leading razor protein`)
        ][,Midpoint:= map_dbl(.x = Sequence,~str_locate(tmp_prot_seq,.x) %>% mean())] 
        
        
        
        annotation <- tmp_pep_col[,.(Sequence,Proteins,`Leading razor protein`, Midpoint)] %>% column_to_rownames("Sequence") %>% 
            set_names(c("PG","Majority","Peptide Position")) %>% 
            mutate(
                # Isoforms = case_when(
                # str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4")~"Q9UPN3-5;Q9UPN3-4",
                # str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4",negate = T)~"Q9UPN3-5",
                # str_detect(Isoforms,"-5",negate = T) & str_detect(Isoforms,"-4")~"Q9UPN3-4",        ),
                `Peptide Position` = rank(`Peptide Position`))
        protein_order <- annotation %>% subset(rownames(.) %in% rownames(tmp_heatmap)) %>% arrange(`Peptide Position`) %>% rownames()
        # tmp <- tmp %>% as.matrix %>% .[names_order,names_order]
        # tmp[is.na(tmp)] <- 0
        pheatmap::pheatmap(tmp_heatmap[protein_order,protein_order],show_rownames = F,show_colnames = F,
                           cluster_rows = F,cluster_cols = F,
                           annotation = annotation, color=myColor, breaks=myBreaks, main = tmp_protein)
    }
    
    # testing <- amap::Dist(tmp,"correlation")
    
    
    
    tmp <- as.data.table( reshape2::melt( as.matrix( tmp  )))
    tmp <- tmp[, .( Peptide1 = as.character(Var1), Peptide2 = as.character(Var2), PCC = value ) ]
    tmp <- tmp[ Peptide1 >= Peptide2 ][!is.nan(PCC)] #%>% distinct()
    #if using leading razor, then we can then look back at a simplified evidence table to see if the peptides were unique or shared in Proteins
    tmp[, Protein := tmp_protein ]
    correlations <- rbind( correlations, tmp)
    setTxtProgressBar(pb, i)
}

tmp_protein <- "Q03001" #unique(med_DT$`Leading Razor Protein`)[i] #"Q9UPN3"# 
tmp <- med_DT[ `Leading Razor Protein` == tmp_protein ]
# Should we log 2 the median Ratio?  now the positive FC might have a larger influence on corr than negative
tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
tmp$Experiment <- NULL
#finding and removing correlations with less than 15 common datapoints
missingness_tmp <- crossprod(is.na(tmp))>10
# Distance <- amap::Dist(t(tmp),method = "correlation",diag = TRUE, upper = TRUE)
# Distance[missingness_tmp == FALSE] <- 1
# Distance[is.na(Distance)] <- 1
# Distance[is.na(Distance)] <- 1
tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")#using suggested bicor
tmp[missingness_tmp == FALSE] <- 0
tmp[is.na(tmp)] <- 0
is_0 <- tmp==0
keep_few_0 <- is_0 %>% matrixStats::rowSums2() %>% set_names(rownames(tmp)) %>% .[.<(ncol(tmp)/(1.3))] %>% names() #only keep peptides with few 0/NA
tmp <- tmp[keep_few_0,keep_few_0]
# res <- dbscan::optics(as.dist(1-tmp), minPts = 5)
# res <- dbscan::extractXi(res, xi = 0.05)
# res$clusters_xi

# hdbscan_peptides <- dbscan::hdbscan(Distance, minPts = 5)
# annotation <- data.frame(Cluster = res$cluster
#                          # hdbscan_peptides$cluster
#                          ,
#                          Seq = colnames(tmp) )%>% column_to_rownames("Seq")
# names_order <- annotation %>% arrange(Cluster) %>% rownames()
tmp_pep_col <- Peptides_cols[`Leading razor protein`==tmp_protein]

annotation <- tmp_pep_col[,.(Sequence,Proteins,Midpoint)] %>% column_to_rownames("Sequence") %>% 
    set_names(c("Isoforms","Peptide Position")) %>% 
    mutate(
        # Isoforms = case_when(
    #     str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4")~"Q9UPN3-5;Q9UPN3-4",
    #     str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4",negate = T)~"Q9UPN3-5",
    #     str_detect(Isoforms,"-5",negate = T) & str_detect(Isoforms,"-4")~"Q9UPN3-4",
    # ),
    `Peptide Position` = rank(`Peptide Position`))
protein_order <- annotation %>% subset(rownames(.) %in% rownames(tmp)) %>% arrange(`Peptide Position`) %>% rownames()

pheatmap::pheatmap(tmp[protein_order,protein_order],show_rownames = F,show_colnames = F,
                   cluster_rows = F,cluster_cols = F,
                   annotation = annotation, color=myColor, breaks=myBreaks, main = tmp_protein)


tmp_pep_col <- Peptides_cols[`Leading razor protein`=="Q9UPN3"]

annotation <- tmp_pep_col[,.(Sequence,Proteins,Midpoint)] %>% column_to_rownames("Sequence") %>% 
    set_names(c("Isoforms","Peptide Position")) %>% 
    mutate(
        Isoforms = case_when(
            str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4")~"Q9UPN3-5;Q9UPN3-4",
            str_detect(Isoforms,"-5") & str_detect(Isoforms,"-4",negate = T)~"Q9UPN3-5",
            str_detect(Isoforms,"-5",negate = T) & str_detect(Isoforms,"-4")~"Q9UPN3-4",
        ),
        `Peptide Position` = rank(`Peptide Position`))



Q9UPN3_intensity <- DT[str_detect( `Leading Razor Protein`,"Q9UPN3")
                       ][, .(median_Intensity = median(Intensity, na.rm = TRUE)), .(Sequence, Proteins, Experiment, `Leading Razor Protein`)
                         ][, median_Intensity:= log10(median_Intensity)] 
Q9UPN3_intensity <-     dcast(Q9UPN3_intensity, Experiment + Proteins ~ Sequence, value.var = "median_Intensity" ) %>%
    melt(id.vars = c("Experiment","Proteins"), measure.vars = colnames(.)[-c(1:2)]) %>% 
    .[annotation %>% rownames_to_column("Sequence"), on = .(variable= Sequence)]
Q9UPN3_intensity <- Q9UPN3_intensity[,.(Isoform_intensity = median(value, na.rm = T)), by = .(Experiment,Isoforms)][
    !is.na(Experiment)
]
Q9UPN3_intensity[,PXD001406 := fifelse(str_detect(Experiment,"1406"),"PXD001406","else")]
# Q9UPN3_intensity$Experiment <- factor(Q9UPN3_intensity$Experiment,
#                                       Q9UPN3_intensity[Isoforms == "Q9UPN3-4"][order(Isoform_intensity),Experiment])
# 
    # .[annotation,on = .(variable= Sequence)] %>%
        # .[1:1000,] %>% 
    ggplot(Q9UPN3_intensity,aes(x  = Experiment, y = Isoform_intensity, group = Isoforms, colour = Isoforms))+
    geom_point(alpha = 0.5)+
    geom_line()+
        theme_bw()+
    # theme(legend.position = "none")+
        facet_wrap("PXD001406",scales = "free_x")+
    ggtitle("Q9UPN3",
            subtitle = "PXD001406	EBV-transformed lymphoblastoid cell lines (LCLs) from 62 Yoruba HapMap individuals compared to a reference individual
	EBV-transformed lymphoblastoid cell lines (LCLs) from 62 Yoruba HapMap individuals compared to a reference individual")



correlations <- correlations[,n:=.N, by = Protein][n>15]



for_diffacto <- med_DT[str_detect(file, "^g3_", negate = T)] %>% 
    dcast( sequences+Proteins~file, value.var = "median_ratio" )
multi_peptide_proteins <- for_diffacto[,.N,by = Proteins][order(N)][N>15,Proteins]
for_diffacto <- for_diffacto[Proteins %chin% multi_peptide_proteins]
write.csv(for_diffacto,"ProteomeHD_example_pep.csv",#na = "",
          row.names = F)
write_tsv(data.frame(X1 = colnames(for_diffacto)[-c(1:2)],
                                   X2 = str_match(colnames(for_diffacto)[-c(1:2)],"(^[:graph:]*?)_")[,2]),"samples_Proteome_HD.tsv",
          col_names = F)
med_DT <- DT[, .(median_ratio = median(`Ratio H/L normalized`, na.rm = TRUE)), .(Sequence, Proteins, Experiment, `Leading Razor Protein`)]
setnames(med_DT, old = names(med_DT), new = c("sequences","Proteins", "file","Protein","median_ratio") )

i <- unique(med_DT$Protein) %>% str_detect("Q9Y490") %>% which
tmp_protein <- unique(med_DT$Protein)[i]
tmp <- med_DT[ Protein == tmp_protein ]
PLEC_loading <- fread(here::here("Datasets","Raw","ProteomeHD.pep_loading.txt"))[str_detect(V1,"Q9Y490"),!c("V1"), with = FALSE]
PLEC_loading %>% setnames(old = names(.),new = c("Sequence","Diffacto_Loading"))
temp <- tmp[PLEC_loading, on = c("sequences" = "Sequence")][,Quant := fifelse(Diffacto_Loading>0.5,T,F)]    
ggplot(temp %>% na.omit()%>% pivot_wider(names_from = file, values_from = median_ratio) %>% pivot_longer(
    cols = colnames(.) %>% str_subset("^g"), names_to = "file", values_to = "median_ratio"
) %>% subset(!between(Diffacto_Loading, 0.3,0.7)),aes(x  = reorder(file,median_ratio), y = median_ratio, group = sequences, colour = Quant))+
    geom_point(alpha = 0.1)+
    geom_line(alpha = 0.8)+
    ggtitle(tmp_protein)

cl <- hdbscan(amap::Dist(rbind(testing,testing),method = "correlation") , minPts = 5)
cl

#Do common peptides assigned to other leading protein groups by razor better correlate with peptides in another leading protein?

tmp_protein <- "Q9UPN3" #unique(med_DT$`Leading Razor Protein`)[i] #"Q9UPN3"# 
razor_proteins<- med_DT[ str_detect(Proteins,tmp_protein),Proteins ] %>% unique() %>% str_split(";") %>% unlist()%>% unique() %>% 
    str_subset(tmp_protein)
tmp <- med_DT[`Leading Razor Protein` %in% razor_proteins]
# Should we log 2 the median Ratio?  now the positive FC might have a larger influence on corr than negative
tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
tmp$Experiment <- NULL
#finding and removing correlations with less than 15 common datapoints
missingness_tmp <- crossprod(is.na(tmp))>10
#Distance <- amap::Dist(t(tmp),method = "correlation",diag = TRUE, upper = TRUE)
#Distance[missingness_tmp == FALSE] <- 1
#Distance[is.na(Distance)] <- 1
# Distance[is.na(Distance)] <- 1
tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")#using suggested bicor
tmp[missingness_tmp == FALSE] <- 0
tmp[is.na(tmp)] <- 0
is_0 <- tmp==0
keep_few_0 <- is_0 %>% matrixStats::rowSums2() %>% set_names(rownames(tmp)) %>% .[.<(ncol(tmp)/3)] %>% names() #only keep peptides with few 0/NA
tmp <- tmp[keep_few_0,keep_few_0]
# res <- dbscan::optics(as.dist(1-tmp), minPts = 5)
# res <- dbscan::extractXi(res, xi = 0.05)
# res$clusters_xi

# hdbscan_peptides <- dbscan::hdbscan(Distance, minPts = 5)
# annotation <- data.frame(Cluster = res$cluster
#                          # hdbscan_peptides$cluster
#                          ,
#                          Seq = colnames(tmp) )%>% column_to_rownames("Seq")
# names_order <- annotation %>% arrange(Cluster) %>% rownames()
 # tmp_pep_col <- Peptides_cols[`Leading razor protein`==tmp_protein]

annotation <- DT[Sequence %chin% colnames(tmp),.(Sequence,Proteins,`Leading Razor Protein`)] %>% distinct() %>% 
    column_to_rownames("Sequence") %>% 
    set_names(c("Protein_group","Assigned_Protein"))
pept_order <- annotation %>% arrange(Protein_group) %>% rownames()
pheatmap::pheatmap(tmp[pept_order,pept_order],
                   show_rownames = F,show_colnames = F,
                   # cluster_rows = F,cluster_cols = F,
                   annotation = annotation, color=myColor, breaks=myBreaks, main = "Peptides shared with MACF1")

tmp_protein <- "Q9UPN3" #unique(med_DT$`Leading Razor Protein`)[i] #"Q9UPN3"# 
razor_proteins<- med_DT[ str_detect(Proteins,tmp_protein),Proteins ] %>% unique() %>% str_split(";") %>% unlist()%>% unique() 
tmp <- med_DT[`Leading Razor Protein` %in% razor_proteins]
# Should we log 2 the median Ratio?  now the positive FC might have a larger influence on corr than negative
tmp <- dcast(tmp, Experiment ~ Sequence, value.var = "median_ratio" )
tmp$Experiment <- NULL
#finding and removing correlations with less than 15 common datapoints
missingness_tmp <- crossprod(is.na(tmp))>10
#Distance <- amap::Dist(t(tmp),method = "correlation",diag = TRUE, upper = TRUE)
#Distance[missingness_tmp == FALSE] <- 1
#Distance[is.na(Distance)] <- 1
# Distance[is.na(Distance)] <- 1
tmp <-  WGCNA::bicor(tmp,use = "pairwise.complete.obs") #cor(tmp, use = "pair")#using suggested bicor
tmp[missingness_tmp == FALSE] <- 0
tmp[is.na(tmp)] <- 0
is_0 <- tmp==0
keep_few_0 <- is_0 %>% matrixStats::rowSums2() %>% set_names(rownames(tmp)) %>% .[.<(ncol(tmp)/2)] %>% names() #only keep peptides with few 0/NA
tmp <- tmp[keep_few_0,keep_few_0]
# res <- dbscan::optics(as.dist(1-tmp), minPts = 5)
# res <- dbscan::extractXi(res, xi = 0.05)
# res$clusters_xi

# hdbscan_peptides <- dbscan::hdbscan(Distance, minPts = 5)
# annotation <- data.frame(Cluster = res$cluster
#                          # hdbscan_peptides$cluster
#                          ,
#                          Seq = colnames(tmp) )%>% column_to_rownames("Seq")
# names_order <- annotation %>% arrange(Cluster) %>% rownames()
# tmp_pep_col <- Peptides_cols[`Leading razor protein`==tmp_protein]

annotation <- DT[Sequence %chin% colnames(tmp),.(Sequence,Proteins,`Leading Razor Protein`)] %>% distinct() %>% 
    column_to_rownames("Sequence") %>% 
    set_names(c("Protein_group","Assigned_Protein"))
pept_order <- annotation %>% arrange(Protein_group) %>% rownames()
pheatmap::pheatmap(tmp[pept_order,pept_order],
                   show_rownames = F,show_colnames = F,
                   # cluster_rows = F,cluster_cols = F,
                   annotation = annotation, color=myColor, breaks=myBreaks, main = "Peptides shared with MACF1")

