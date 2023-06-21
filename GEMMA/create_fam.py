import pandas as pd
import sys
fam_path = sys.argv[1]
pheno_path = sys.argv[2]
fam = pd.read_csv(fam_path, sep=" ", encoding='latin1', header=None)
fam.iloc[:,5] = 1 #remplace the last column of the maf file (with only missing pheno by only present pheno
#so the relatdeness matrice cpmputaion will happen using all th eindividuals included

pheno = pd.read_csv(pheno_path, sep= "\t", encoding='latin1', header=None)
pheno.drop(labels=1, axis=1, inplace=True) #remove the famid (duplicate of the ID)

merged_df = fam.merge(pheno, how= "left", on=0)
merged_df.to_csv(fam_path, sep=' ', header=None, index=None, na_rep="NA")