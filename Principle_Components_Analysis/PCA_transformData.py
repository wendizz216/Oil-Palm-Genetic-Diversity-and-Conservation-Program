# !/usr/bin/env python

#Python code is for categoring snp/genotype as minor, major or het for each snp position for all genotypes
#This will be used for Multiple regression where the allele frequency of each snp position will be the response
#variable.
#Numerical indicator of genotypic state: 0 = homo minor allele genotype; 0.5 = het; 1=homo major allele genotype
#maf > 0.01; no more than 5% missing genotype call rate for freq estimates


import os
import numpy as np
import scipy as sp
import itertools
from itertools import izip

frq = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/PCA/Eguineensis.VF.SNP.95.frq",dtype="str", delimiter=("\t"),skiprows=1)


positions=[]
freq_1 = []
freq_2 = []

for pos in frq:
	#print pos
	freq_1.append(pos[4].split(":"))
	freq_2.append(pos[5].split(":"))
	positions.append(pos[0:4])


freq = np.hstack((np.asarray(positions),np.asarray(freq_1), np.asarray(freq_2)))

##Filtered snp positions from frequency file 

geno=np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/Eguineensis.VF.SNP.95.vcf.tab",
                dtype="str",delimiter="\t",skiprows=1)

samples=np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/sample_names.txt",
                dtype="str",delimiter="\t")

amer_sur_g = np.loadtxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/PCA/amer_sur_g_samples.txt", dtype="str",delimiter="\t")

#print samples

pop_idx = [i for i, item in enumerate(list(samples)) if item in list(samples)]

print pop_idx, len(pop_idx)

geno_header=samples[pop_idx]
#print geno_header
geno = geno[:,pop_idx]

freq1=np.array([x[0] + "_" + x[1] for x in freq])
geno1=np.array([x[0] + "_" + x[1] for x in geno])

#print freq1, geno1


indices=np.nonzero(np.in1d(geno1,freq1)) #gives you indices of the first array --> geno1

#print indices

geno_frq = geno[indices] #positions that intersect with the frequency file "frq"

##This code takes the first allele (A1) freq and determines major and minor and hets relative to A1 allele

allele_id=[]
for pos in freq:
    if float(pos[5]) < float(pos[7]):
        allele_id.append("minor")
    	if float(pos[5]) == float(pos[7]):
		allele_id.append("moderate")
    else:
        allele_id.append("major")

#print allele_id

allele_type=[x.split(",") for x in allele_id] #make list into an 1d array to append columns
#print np.asarray(allele_type)

frq_snp=np.hstack((freq,np.asarray(allele_type)))

#print len(freq) #This is the array with major and minor allele indicator column
#print freq_snp


#These are a lists of genotypes that indicate major, minor alleles
major=[]
minor=[]
het1=[]
het2=[]

for pos in frq_snp:
    if pos[8] == "major":
        major.append([pos[4]+"/"+pos[4]])
        minor.append([pos[6]+"/"+pos[6]])
        het1.append([pos[4]+"/"+pos[6]])
        het2.append([pos[6]+"/"+pos[4]])
    
    elif pos[8] == "minor":
        major.append([pos[6]+"/"+pos[6]])
        minor.append([pos[4]+"/"+pos[4]])
        het1.append([pos[4]+"/"+pos[6]])
        het2.append([pos[6]+"/"+pos[4]])

    elif pos[8] == "moderate":
        major.append([pos[4]+"/"+pos[4]])
        minor.append([pos[6]+"/"+pos[6]])
        het1.append([pos[4]+"/"+pos[6]])
        het2.append([pos[6]+"/"+pos[4]])
        
#print len(major),len(minor)
#print major, minor


#Column header:
#["chrom","pos","N_Alleles","N_chrom","A1","freq_A1","A2","freq_A2","Freq_cat","major","minor","het1","het2"]
      
genetic_state=np.hstack((freq,np.asarray(major),np.asarray(minor),np.asarray(het1),np.asarray(het2)))

#print genetic_state


genotype_state_file=[]

#print geno_frq
#print genetic_state

for pos, cat in itertools.izip(geno_frq, genetic_state):        
        for n, i in enumerate(pos):
            #print n, i

            if i == cat[8]:
                pos[n]=1 # Major allele = 1

            elif i == cat[9]:
                pos[n]=0 # Minor allele = 0

            elif i == "./.":
                pos[n]="NA" #Missing genotypes
        
            elif i == cat[10]:
                pos[n]=0.5 # Het = 2

            elif i == cat[11]:
                pos[n]=0.5 #repeat to ensure all possible hetero genotypes are identified
	genotype_state_file.append(pos)


final=np.vstack([geno_header,genotype_state_file])
        
print final, len(geno_header)

s = 189*"%s,%s,"+"%s,%s"

np.savetxt("/home/cmb-07/sn1/wvu/Oil_Palm/analysis_April2016/EGuineensis/PCA/Eguineensis.VF.SNP.95.Major_Minor_Allele_all.txt", final, fmt=s,delimiter=",")

