# Archaic introgression for candidate regions for selection by hmmix (Skov et al., 2018)


there are two folders, **run_hmmix** and **process_res**: <br>

A) **run_hmmix** contains the following scripts and implements the clear workflow that could be found here <br>
https://github.com/LauritsSkov/Introgression-detection <br>

1. processing_gvcf.sh # bcftools commands that were used for processing the (g)vcf files for the hmmix <br>
2. Skov_pipeline.preparing.sh # preparing specific input for the hmmix <br>
3. unite_chr_per_sample.sh # combining input for the hmmix model <br>
4. Skov_pipeline.hmm.sh # training and decoding by the hmmix model <br>
5. laucnh_Skov_pipeline.sh # launching 2-4 bullet points in sbatch mode for the regions of interest <br>

B) **process_res** contains the following scripts for proccessing the hmmix output <br>

Introgression_processing.workflow.R # includes the workflow for finding the top frequent haplotype per candidate region for selection <br>
Pairwise_comparison.chr1_GBP.R # includes the counting of the pairwise distance between archaics in GBP region on chr1 <br>
