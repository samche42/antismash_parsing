#Usage python domain_extractor.py /path/to/antismash/gbk_files

from Bio import SeqIO
import pandas as pd
import sys
import os
import numpy as np

input_directory = sys.argv[1]

gbk_files = [file for file in os.listdir(input_directory) if file.endswith(".gbk")]

domain_types = ["ACP","ACP_beta","ACPS","Aminotran_1_2","Aminotran_3","Aminotran_4","Aminotran_5","AMP-binding","A-OX","B","CAL_domain","Cglyc","cMT","Condensation","Condensation_DCL","Condensation_Dual","Condensation_LCL","Condensation_Starter","ECH","Ene_KS","Epimerization domain","F","FkbH","GNAT","Hal","Heterocyclization","Hyb_KS","Itr_KS","Mod_KS","NAD_binding_4","nMT","NRPS-COM_Cterm","NRPS-COM_Nterm","oMT","PCP","PKS_AT","PKS_DH","PKS_DH2","PKS_DHt","PKS_Docking_Cterm","PKS_Docking_Nterm","PKS_ER","PKS_KR","PKS_KS","PKS_PP","Polyketide_cyc2","Polyketide_cyc","PS","SAT","TD","Thioesterase","TIGR01720","TIGR02353","Tra_KS","Trans-AT_docking","X"]

for domain in domain_types:
	output_string = ''
	dom = str(domain)
	for gbk_file in gbk_files:
		count = 0
		BGC = SeqIO.read(gbk_file,'genbank')
		file_name = str(gbk_file)
		molecule = file_name.replace('.gbk','')
		for feat in BGC.features:
			if feat.type == 'aSDomain':
				if feat.qualifiers['aSDomain'][0] == dom:
					count +=1
					descrip = feat.qualifiers['aSDomain'][0]
					gene_prot = feat.qualifiers['translation'][0]
					output_string = output_string+">"+molecule+"_"+descrip+str(count)+"\n"+gene_prot+"\n"
	if len(output_string) > 0:
		outputfile = str(domain)+"_seqs.faa"
		with open(input_directory+"/"+outputfile, "w") as text_file:
			text_file.write(output_string)
	else:
		pass
