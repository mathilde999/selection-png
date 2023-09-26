setwd("~/selection-png/Archaic_introgression/process_res/")

library(data.table)
library(dplyr)

## reading all concatenated hmmix results,

## in format 
##[1] "chrom"               "start"               "end"                 "haplo_Skov"          "length"              "state"              
##[7] "mean_prob"           "snps"                "admixpopvariants"    "AltaiNeandertal"     "Chagyrskaya-Phalanx" "Denisova"           
##[13] "Vindija33.19"    

db_Skov <- fread("hmmix.introgression_for.regions_of_selection.raw.tsv.gz")

##### PART 1: filtering the raw results

head(db_Skov)
table(db_Skov$state)
db_Skov <- db_Skov[db_Skov$state == "Archaic",] # step_1: keeping only Archaic
db_Skov$chr_st_end_Skov <- paste(db_Skov$chrom,db_Skov$start,db_Skov$end, sep = ".")
table(db_Skov$state)
db_Skov <- db_Skov[db_Skov$end-db_Skov$start != 0,] # step_2: removing artifacts with explicit segment length == 0


db_Skov2 <- db_Skov[db_Skov$mean_prob >= 0.8 & db_Skov$admixpopvariants >= 1,] #step_3: filtering by mean_prob and archaic variants
str(db_Skov2)

## in this set-up, analysis was launched and concatenated for all samples in the dataset, 
## here we're keeping only HL and LL PNG populations
smlp <- fread("High_Low.ID", header = F)$V1 
db_Skov2$sample <-gsub('.{5}$', '', db_Skov2$haplo)

db_Skov2 <- db_Skov2[db_Skov2$sample %in% smlp,]


## in this set-up, we considered more vise to launch it to the whole chromosome of the interest, 
## so here we need to subset the candidate regions

write.table(db_Skov2[order(db_Skov2$chrom,db_Skov2$start,db_Skov2$end),1:3],"hmmix.filtered.bed", col.names = F, row.names = F, sep = "\t", quote = F)

beds <- fread("regions_and_SNP_of_interest.csv") # process candidates
beds$CHR <- as.numeric(gsub("^chr","",beds$CHR))
str(beds)

write.table(beds[order(beds$CHR,beds$START,beds$END),1:3],"candidates.44.bed", col.names = F, row.names = F, sep = "\t", quote = F)

# eventhough it could be launched via R, the person in charge still prefers cmd tool
system("bedtools intersect -a candidates.44.bed -b hmmix.filtered.bed -wao > overlap.bed") 


# preparing for merging and merge itself of 3 tables
Skov2_bedtools <- fread("overlap.bed") %>% unique()
colnames(Skov2_bedtools) <- c("CHR","START","END","chrom","start","end","ovelap")
Skov2_bedtools$length_Skov <- Skov2_bedtools$end-Skov2_bedtools$start
Skov2_bedtools$ratio <- Skov2_bedtools$ovelap/Skov2_bedtools$length_Skov

Skov2_bedtools$chrom <- as.numeric(Skov2_bedtools$chrom)
Skov2_bedtools[is.na(Skov2_bedtools$chrom),]
Skov2_bedtools <- Skov2_bedtools[!is.na(Skov2_bedtools$chrom),] 

str(db_Skov2)

res <- merge(Skov2_bedtools,db_Skov2) %>% merge(.,beds, by = c("CHR","START","END")) %>% unique()

quantile(res$ratio)

res <- res[res$ratio == 1,] # step_4: keeping only those that fully fall into the candidate region
res[,1:3] %>% unique() # which region out of 44 got some intrgoressed haplotypes



##### PART 2: annotating filtered individual-level results 

res$plain_origin <- NA

res$plain_origin[res$Denisova > res$AltaiNeandertal & 
                   res$Denisova > res$`Chagyrskaya-Phalanx` &
                   res$Denisova > res$Vindija33.19] <- "Denisova" 


res$plain_origin[res$Denisova < res$AltaiNeandertal & 
                   res$Denisova < res$`Chagyrskaya-Phalanx` &
                   res$Denisova < res$Vindija33.19] <- "Neandertal" 


res$plain_origin[is.na(res$plain_origin)] <- "ambiguous"
res$plain_origin %>% table()

# straight comparison of the manual annotating vs hmmix one. it's noisy since it's still on individual-level haplotypes, 
# and due to the complex structure of PNG regarding the introgression
res %>% select(Introgression,plain_origin) %>% table()  


high <- fread("High.ID", header = F)$V1
low <- fread("Low.ID", header = F)$V1

res$pop[res$sample %in% high] <- "HL"
res$pop[res$sample %in% low] <- "LL"
sum(is.na(res$pop))
res$pop %>% table()

write.table(res,"hmmix_segements.per_sample.per_regions.filtered_conservative.tsv", quote = F, row.names = F, sep = "\t")



##### PART 3: collapsing individual-level results into regions of introgression


## gathering from observation data all aSNPs; that's hard to share this data, so we share the resulting table only
## nothing important here, rather time-consuming and not optimized subsetting
res$chrom <- gsub("^","chr",res$chrom)

#df_tot <- data.frame()
#for (i in 1:nrow(res)){
#  name <- res$sample[i]
#  df <- fread(paste0("../run_hmmix/out_gvcf_filt_phased/observ/obs.",name,".txt"))
#  df <- df[df$chrom == res$chrom[i] & df$pos > res$start[i] & df$pos < res$end[i],]
#  df_tot <- rbind(df_tot,df)
#}
#db_colapsed <- df_tot[order(df_tot$pos),1:2] %>% unique()
#write.table(db_colapsed[order(db_colapsed$chrom,db_colapsed$pos),],"snps.res.all.tsv", row.names = F, col.names = F, sep = "\t", quote = F)


## counting R2 between the overmentioned aSNPs, and collapsing them into raw regions of introgression based on R2 >= 0.5 cut-off
## done in cmd via bcftools + vcftools

#for i in {1..4} {6..10} {12..14} {17..19} 22; do
#bcftools view -R snps.res.all.tsv -S High_Low.ID ../run_hmmix/vcfs_all_positions/phasing/chr${i}_phased.YRI_PAP_gen_vcf.bcf -Ob |\
#vcftools --bcf - --hap-r2 --min-r2 0.05 --out auxiliary_r2/chr$i && bgzip -c auxiliary_r2/chr$i.hap.ld > auxiliary_r2/chr$i.hap.ld.gz && rm auxiliary_r2/chr$i.hap.ld; done

#cd auxiliary_r2; mdkir res

#for i in {1..4} {6..10} {12..14} {17..19} 22; do Rscript simplets.collapse.R $i; done
#for i in res/chr*haplotypes.R2_0_5.tsv; do tail -n+2 $i >> haplotypes.R2_0_5.tsv; done
#for i in res/chr*.aSNPs.R2_0_5.tsv; do tail -n+2 $i >> aSNPs.R2_0_5.tsv; done


## processing the raw results

db <- fread("auxiliary_r2/haplotypes.R2_0_5.tsv")

## to decrease the number of spurious regions, 
## we restricted our analysis to inferred regions of introgression longer than 15kb and having more than or equal to 5 aSNPs
db2 <- db[db$N_aSNP >= 5 & db$length > 15000,]

df <- fread("auxiliary_r2/aSNPs.R2_0_5.tsv")

db2 <- merge(df,db2)
db2$chr_st_end %>% unique() %>% length()
colnames(db2)[2:3] <- c("chrom","pos")

db2_gr <- db2 %>% group_by(chr_st_end) %>% 
  summarise(total_in_reg = n())

df_res <- data.frame()

## the complex function to find overlap between aSNP linked by R2 >= 0.5 (regions of introgression)
## and the aSNPs in individual introgressed haplotypes

#gettbl <- function(i){
#  name <- res$sample[i]
#  df <- fread(paste0("../run_hmmix/out_gvcf_filt_phased/observ/obs.",name,".txt"))
#  df <- df[df$chrom == res$chrom[i] & df$pos > res$start[i] & df$pos < res$end[i],]
#  df_gr <- nrow(df) # total not splits (need all)
  
#  df_db <- merge(df,db2, by = c("chrom","pos"))
#  df_db_gr <- nrow(df_db)
  
#  df_db <- df_db %>% group_by(chr_st_end) %>% 
#    summarise(total_overlap=n())
  
#  df_db <- merge(df_db,db2_gr)
  
#  met <- data.frame(sample = res$sample[i],
#                    haplo_Skov = res$haplo_Skov[i],
#                    chrom = res$chrom[i],
#                    start = res$start[i],
#                    end = res$end[i],
#                    count_across_all_reg = df_gr,
#                    ovrlp_across_all_reg = df_db_gr)
#  if (nrow(df_db) == 0){
#    ttl <- cbind(met,data.frame(chr_st_end = NA,
#                                total_overlap = NA,
#                                total_in_reg = NA))
#  } else {
#    ttl <- cbind(met,df_db)
#  }
#  
#  
#  return(ttl)
#}

## and the function applying

#new_res <- lapply(1:nrow(res), function(i) gettbl(i)) %>% bind_rows()
#write.table(new_res,"hmmix.segements_per_sample.grouping.tsv", quote = F, row.names = F, sep = "\t")

new_res <- fread("hmmix.segements_per_sample.grouping.tsv")
new_res2 <- new_res[!is.na(new_res$chr_st_end),]

new_res2$rat <- new_res2$total_overlap/new_res2$total_in_reg

quantile(new_res2$rat,seq(0,1,0.01))

## since some individuals could carry recombinant haplotypes
## we apply a cut-off to count a haplotype as present if it 
## i) its aSNPs overlap with aSNPs in the region of introgression with 0.3 ratio 
## ii) raw number of aSNPs from overlap larger than 5, in order to filter stronger short haplotypes
## the cut-off helps to decrease false positive matches
new_res3 <- new_res2[new_res2$rat > 0.3 & new_res2$total_overlap >= 5,c(1,2,3,4,5,8)] %>% unique()


## merging tables 
res2 <- merge(res,new_res3, by = c("sample","haplo_Skov","chrom","start","end"))
res2$chr_st_end %>% unique() %>% length()

res2 <- res2[order(res2$chr_st_end),]
write.table(res2,"hmmix.uncollapsed_filtered.with_grouping.tsv", quote = F, row.names = F, sep = "\t")


##### PART 4: counting haplotypes' frequencies in regions of introgression

## counting frequencies and ancestries in the region per population
sums_res <- res2 %>% group_by(chr_st_end,CHR,START,END,pop) %>% 
  summarise(Introgression_vec=paste(unique(sort(Introgression)),collapse = ","),
            plain_origin_vec=paste(unique(sort(plain_origin)),collapse = ","),
            haplo_count=n())

## counting ancestry frequencies, ancestry counts
sums_tmp <- res2 %>% group_by(chr_st_end,CHR,START,END,plain_origin) %>% 
  summarise(orig_counts=n()) %>% filter(orig_counts == max(orig_counts))

## hotfix for two rare haplotypes of "ambiguous" category
sums_tmp[duplicated(sums_tmp[,1:4]) | duplicated(sums_tmp[,1:4], fromLast=TRUE),]
sums_tmp$plain_origin[duplicated(sums_tmp[,1:4]) | duplicated(sums_tmp[,1:4], fromLast=TRUE)] <- "ambiguous"
sums_tmp <- sums_tmp %>% unique()

## counting ancestry frequencies, total counts
sums_tmp2 <- res2 %>% group_by(chr_st_end,CHR,START,END) %>% 
  summarise(total_counts=n())

## counting ancestry frequencies, merging and counting frequency
sums_tmp <- merge(sums_tmp,sums_tmp2)
sums_tmp$category_ratio <- sums_tmp$orig_counts/sums_tmp$total_counts
sums_res_tmp <- merge(sums_res,sums_tmp[,c(1:5,8)])


## choosing the most prevalent haplotype per population per candidate region of selection
sums_max <- sums_res_tmp %>% group_by(CHR,START,END,pop) %>% 
  filter(haplo_count == max(haplo_count))

write.table(sums_max,"hmmix.top_freq_haplotypes.per_candidate_region.tsv", quote = F, row.names = F, sep = "\t")

## uncollapsed data for top_freq_haplotypes

uncol_top <- res2[res2$chr_st_end %in% sums_max$chr_st_end,]
uncol_top <- uncol_top[order(uncol_top$chr_st_end),]

write.table(uncol_top,"hmmix.top_freq_haplotypes.per_candidate_region.uncollapsed.tsv", quote = F, row.names = F, sep = "\t")


uncol_2 <- res2[res2$chr_st_end %in% c("chr1.89054418.89202534","chr12.58451248.58568114"),]
uncol_2 <- uncol_2[order(uncol_2$chr_st_end),]

write.table(uncol_2,"hmmix.two_haplos.uncollapsed.tsv", quote = F, row.names = F, sep = "\t")
