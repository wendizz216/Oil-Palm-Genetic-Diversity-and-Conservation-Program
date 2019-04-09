# !/usr/bin/env python

#Sampling individuals across different sample sizes to identify the least number of genotypes that capture the maximum genetic diversity. here we use number of segregating sites to assess genetic diversity. A genetic diversity index will determine the proportion of genetic diversity captured across different samples sizes to identify the minimal sample size. genetic diversity index == number of segregating sites at N sample size/number of segregating sites all samples

import numpy as np
import pandas as pd
import random
import os, sys
import fileinput
import glob
import itertools


tab = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/SOL_OILPALM.Eguineensis9.1.VF.SNP.DP8.90.vcf.maxAllele2.tab",dtype='str',delimiter="\t",skiprows=1)

samples = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/QC/sample_names.txt",dtype="str",delimiter="\t")
hcpc_tp_bra_sur = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/samples/hcpc_tp_bra_sur.txt",dtype="str",delimiter="\t")

hcpc = [s for s in hcpc_tp_bra_sur if any(xs in s for xs in ['HND','CRI','PAN','COL'])]


#STEP ONE
#############Create genotype/individual list from random sampling across different sample sizes 100 times
x1=5
x2=150
sample_size = range(x1,x2+5)

for N in sample_size:
        print N
        for ID in range(1,101): #number of iterations/sampling
                print ID

		random.shuffle(hcpc) 

		list_ = []
                list_size = N
		
		for i in range(0, list_size):
			list_.append(hcpc[i])

                output_directory = "/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/diversity"

                list_file='%s/sampling_N/hcpc_N%s_%s' %(output_directory,N,ID) #CHANGE

                out = open(list_file, "w")
                for i in range(0,len(list_)):
                        out.write('%s\n' %(list_[i]))

                out.close()

#########################################################################################################

#STEP TWO
#############Calculate the proportion of genomic variant sites for each sample size (number of genomic variant sites for sample size X/total genomic variant site in all samples, N=150).

#retreive genotype files from sampling directory
source_files = "/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/diversity/sampling_N/*"

missing_rate = 0.049
df = []
for input_file in glob.glob(source_files):
	print input_file
	file_name = input_file.split("sampling_N/")[1]
	group = file_name.split('_')[0]
	N_size = file_name.split('_')[1].replace('N','')
	iteration = file_name.split('_')[2]

	sampled_geno = np.loadtxt("%s" %input_file, dtype="str",delimiter="\t")
	#print sampled_geno
	geno_idx = [i for i, item in enumerate(list(samples)) if item in list(sampled_geno)]
	geno_snp = tab[:,geno_idx]

	n=0
	m=0
	k=0
	l=0
	for snp in geno_snp:
        	string = "".join(snp)
       		new_str = string.replace("/","")
        	tot_alleles = len(new_str)
        	num_miss = new_str.count('.')
        	prop_miss = num_miss/float(tot_alleles)
        	#print snp, prop_miss

		
        	if prop_miss < missing_rate: #set missing genotype threshold
                	cleaned_snp = [x for x in new_str if "." not in x]
                	counter={}
                	for i in cleaned_snp: counter[i] = counter.get(i,0) + 1
                	frq = sorted([ (freq,allele) for allele, freq in counter.items() ], reverse=True)[:3]
                	#print frq
                	n+=1
			if len(frq) > 1:
				#print frq
				m+=1
				if frq[1][0] == 1: #singleton
					#print frq
					k+=1
				if frq[1][0]/float(frq[0][0]+frq[1][0]) >= 0.05: #maf > 0.05
					#print frq
					l+=1


	#print "num_snp_cov100:",n,"num_seg_cov100:", m, "num_seg_singe:",k, "maf>0.05:",l
	prop_seg = float(m)/n*100
	prop_single = float(k)/n*100
	prop_maf05 = float(l)/n*100
	
	df.append([file_name,group,N_size,iteration,N,missing_rate,prop_seg,prop_single,prop_maf05])

print np.asarray(df)

print "num_snp_all:",N,"num_seg_all:", M, "num_seg_singe:",K, "maf>0.05:",L

header='file_name\tgroup\tN_size\titeration\ttot_loci\tmissing_rate\tprop_SegSNPs\tprop_SegSingle\tprop_SegMaf05'
s='%s\t%s\t'*4+'%s'

np.savetxt('/home/cmb-07/sn1/wvu/Oil_Palm/analysis_Oct2017/diversity/num_SegSites_Sampling_miss05.txt',np.asarray(df),fmt=s,header=header)








