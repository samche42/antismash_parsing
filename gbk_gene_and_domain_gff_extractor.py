#Use: python gbk_gene_and_domain_gff_extractor.py /path/to/antismash/gbk_files

from Bio import SeqIO
import pandas as pd
import sys
import os
import numpy as np

input_directory = sys.argv[1]

gene_df = pd.DataFrame(columns=['molecule', 'gene', 'start', 'end','strand','direction'])
subgene_df = pd.DataFrame(columns=['gene', 'subgene', 'from', 'to'])

gbk_files = [file for file in os.listdir(input_directory) if file.endswith(".gbk")]

for gbk_file in gbk_files:
        BGC = SeqIO.read(gbk_file,'genbank')
        file_name = str(gbk_file)
        molecule = file_name.replace('.gbk','')
        for feat in BGC.features:
                if feat.type == 'CDS':
                        gene = feat.qualifiers['locus_tag'][0]
                        location_raw = str(feat.location)
                        location = location_raw.replace('[','').replace(']','').replace('(',':').replace(')','')
                        start, end, direc = location.split(':')
                        if direc =='+':
                                strand = 'forward'
                                direction = '1'
                        else:
                                strand = 'reverse'
                                direction= '-1'
                        gene_df = gene_df.append({'molecule':molecule, 'gene':gene, 'start':start, 'end':end,'strand':strand, 'direction':direction},ignore_index=True)
                elif feat.type == 'aSDomain':
                        gene = feat.qualifiers['locus_tag'][0]
                        descrip = feat.qualifiers['aSDomain'][0]
                        location_raw = str(feat.location)
                        location = location_raw.replace('[','').replace(']','').replace('(',':').replace(')','')
                        start, end, direc = location.split(':')
                        if direc =='+':
                                strand = 'forward'
                        else:
                                strand = 'reverse'
                        subgene_df = subgene_df.append({'gene':gene, 'subgene':descrip, 'sub_start':start, 'sub_end':end},ignore_index=True)
                        
gene_df.to_csv(input_directory+'/genes.tab', sep='\t', index=False)
subgene_table = pd.merge(subgene_df,gene_df,on='gene', how='inner')
subgene_table = subgene_table[['molecule', 'gene', 'start', 'end','strand','direction','subgene', 'sub_start', 'sub_end']]
subgene_table.to_csv(input_directory+'/domains.tab', sep='\t', index=False)
