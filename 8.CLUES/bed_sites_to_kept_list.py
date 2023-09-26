import pandas as pd
import argparse
import glob

parser = argparse.ArgumentParser(description='')
parser.add_argument("-b", "--bed", help="bed file of the regions of interest CHR\tstart\tend")
parser.add_argument("-t", "--tree", help="extracted tree directory")
parser.add_argument("-o", "--out", help="two column file of the sites in the regions of interest present in relate output CHR\tsite")
args = parser.parse_args()

"""
Get input ready
"""
dir_tree = args.tree
bed = open(args.bed)
"""
create the dictionary of sites list per chr (included in the regions of interest 
First loop save the regions of interest in a dictionary (regions_dic)
Second loop: for each regions in the dictionary, will store the sites that were kept in the mut file of the target pop (chr_to_snp)
"""

regions_dic = {}
CHR=[]
for line in bed:
    chr, start, end = line.split()
    chr = chr.replace("chr","")
    if chr not in CHR:
        CHR.append(str(chr))
        regions_dic[chr] = []
    regions_dic[chr].append(int(end))

kept_list = []
not_kept_list = []
for chr in CHR:
    mut_file = [i for i in glob.glob(f"{dir_tree}/*chr{chr}.mut.gz")][0]
    df = pd.read_csv(mut_file, sep=';', usecols=["pos_of_snp"])
    mut_list = list(df['pos_of_snp'])
    for site in regions_dic[chr]:
        if site in mut_list:
            line = " ".join([str(chr), str(site)])
            kept_list.append(line)
        elif site not in mut_list:
            line = " ".join([str(chr), str(site)])
            not_kept_list.append(line)

output_kept = args.out+"_kept.out"
output_not_kept = args.out+"_not_kept.out"

with open(output_kept, 'w') as fp:
    for item in kept_list:
        fp.write("%s\n" % item)

with open(output_not_kept, 'w') as fp:
    for item in not_kept_list:
        fp.write("%s\n" % item)
