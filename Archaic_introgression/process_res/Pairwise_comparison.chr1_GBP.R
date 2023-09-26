setwd("~/selection-png/Archaic_introgression/process_res/")

library(data.table)
library(dplyr)

# candidate region
#bcftools view -r chr1:88800562-89326878 ../archaic_genoms/four_hg38/chr1.4_genomes.hg38.fixed_ref.masked.bcf | bcftools filter -e 'AC == "."' | bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\t%AC\t%AN\t[%GT\t]\n' | bgzip -c > archaics.GBP.all_positions.tsv.gz


## reading data and subsetting
dbarc <- fread("archaics.GBP.all_positions.tsv.gz")
(dbarc <- dbarc[,-11])
colnames(dbarc) <- c("chrom","pos","REF","ALT_arch","AC","AN",fread("Arch.smpl.vcf_order",header = F)$V1)

# haplotype region: chr1.89054418.89202534
db_reg <- dbarc[dbarc$chrom == "chr1" & dbarc$pos >= 89054418 & dbarc$pos <= 89202534,]


## small rearrangments

colnames(db_reg)[7:10]
colnames(db_reg)[7:10] <- c("AltaiNeandertal_old","Chagyrskaya_Phalanx_old","Denisova_old","Vindija33.19_old")

## removing any ./. positions
db_reg <- db_reg[apply(db_reg[,7:10], 1, function(x)(!("./." %in% x))),]

## pairwise comparison table

pairs_all <- list(c("Denisova_old","AltaiNeandertal_old"),
                  c("Denisova_old","Vindija33.19_old"),
                  c("Denisova_old","Chagyrskaya_Phalanx_old"),
                  c("AltaiNeandertal_old","Vindija33.19_old"),
                  c("AltaiNeandertal_old","Chagyrskaya_Phalanx_old"),
                  c("Vindija33.19_old","Chagyrskaya_Phalanx_old"))

## the comparison function
get_pair_comp <- function(i,pairs){
  first_p <- pairs[[i]][1]
  second_p <- pairs[[i]][2]
  cols_ofint <- db_reg %>% select(all_of(pairs[[i]]))
  #cols_ofint <- cols_ofint[cols_ofint[[first_p]] != cols_ofint[[second_p]],]
  
  cols_ofint <- cols_ofint[!(cols_ofint[[first_p]] == "./." | 
                               cols_ofint[[second_p]] == "./."),]
  
  
  cols_ofint <- data.frame(lapply(cols_ofint, function(x) plyr::revalue(x, c("1/1" = "2","0/1"="1","1/0"="1","0/0" = "0"))))
  cols_ofint[[first_p]] <- as.numeric(as.character(cols_ofint[[first_p]]))
  cols_ofint[[second_p]] <- as.numeric(as.character(cols_ofint[[second_p]]))
  nums <- sum(abs(cols_ofint[[first_p]] - cols_ofint[[second_p]]))/(nrow(cols_ofint)*2)
  nums_total <- sum(abs(cols_ofint[[first_p]] - cols_ofint[[second_p]]))
  
  return(data.frame(arch_1 = first_p,
                    arch_2 = second_p,
                    absolut_difference = nums_total,
                    leng_comp = nrow(cols_ofint),
                    ratio_difference = nums))
  
}


## launch the function + renaming

res2 <- lapply(1:length(pairs_all), function(x) get_pair_comp(x,pairs_all)) %>% bind_rows()
res2$arch_1 <- gsub("_old","",res2$arch_1)
res2$arch_2 <- gsub("_old","",res2$arch_2)

write.table(res2,"pairw_comparison.chr1.89054418.89202534.tsv", sep = "\t", quote = F, row.names = F)
