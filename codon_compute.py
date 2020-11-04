#!/usr/bin/env python3

import os, gzip, itertools

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))
# calculate the total number of genes in Salmonella
gene_num = 0
with gzip.open(file1,"rt") as fh1:
    for line in fh1:
        linestrip = line.strip().split('\t')
        if linestrip[0][0] == ">":
            gene_num += 1
        else:
            continue
print("The total number of genes in Salmonella is {}".format(gene_num))
# calculate the total number of genes in Mycobacterium
gene_num2 = 0
with gzip.open(file2,"rt") as fh2:
    for line in fh2:
        linestrip = line.strip().split('\t')
        if linestrip[0][0] == ">":
            gene_num2 += 1
        else:
            continue
print("The total number of genes in Mycobacterium is {}".format(gene_num2))
# calculate the total length of gens in Salmonella
with gzip.open(file1,"rt") as fh1:
    seqs = dict(aspairs(fh1))
    sum_length1 = 0
    for i in seqs:
        sum_length1 = sum_length1 + len(seqs[i])
    print("The total length of genes in Salmonella fasta file is {}".format(sum_length1))
# calculate the total length of gens in Mycobacterium
with gzip.open(file2,"rt") as fh2:
    seqs = dict(aspairs(fh2))
    sum_length2 = 0
    for i in seqs:
        sum_length2 = sum_length2 + len(seqs[i])
    print("The total length of genes in Mycobacterium fasta file is {}".format(sum_length2))
# calculate the G+C percentage of Salmonella
with gzip.open(file1,"rt") as fh1:
    seqs = dict(aspairs(fh1))
    sequence1 = ""
    for i in seqs:
        sequence1 += seqs[i]
#    print(sequence1)
num_CG = 0
for n in sequence1:
    if n == 'C':
        num_CG += 1
    elif n == 'G':
        num_CG += 1
print("The C+G percentage in Salmonella fasta file is {:.1f}%".format(100*num_CG/sum_length1))
# calculate the G+C percentage of Mycobacterium
with gzip.open(file2,"rt") as fh2:
    seqs = dict(aspairs(fh2))
    sequence2 = ""
    for i in seqs:
        sequence2 += seqs[i]
#    print(sequence2)
num_CG = 0
for n in sequence2:
    if n == 'C':
        num_CG += 1
    elif n == 'G':
        num_CG += 1
print("The C+G percentage in Mycobacterium fasta file is {:.1f}%".format(100*num_CG/sum_length2))
# calculate the total number of condons in each genome of Salmonella
codons1 = (sequence1[n:n+3] for n in range(0,len(sequence1),3))
dict_codons1 = {}
for codon in codons1:
    if codon in dict_codons1: 
        dict_codons1[codon] += 1
    else:
        dict_codons1[codon] = 1
print("Each condons number in Salmonella are shown below: \n{}".format(dict_codons1))
sum_num1 = 0
for i in dict_codons1:
    sum_num1 = sum_num1 + dict_codons1[i]
print("The total number of codons in Salmonella fasta file is {}".format(sum_num1))
# calculate the total number of condons in each genome of Mycobacterium
codons2 = (sequence2[n:n+3] for n in range(0,len(sequence2),3))
dict_codons2 = {}
for codon in codons2:
    if codon in dict_codons2: 
        dict_codons2[codon] += 1
    else:
        dict_codons2[codon] = 1
print("Each condons number in Mycobacterium are shown below: \n{}".format(dict_codons2))
sum_num2 = 0
for i in dict_codons2:
    sum_num2 = sum_num2 + dict_codons2[i]
print("The total number of codons in Mycobacterium fasta file is {}".format(sum_num2))
# Print out table with three columns: Codon, Frequency in Sp1, Frequency in Sp2
d1 = dict_codons1
d2 = dict_codons2
for key1 in d1:
    d1[key1] *= 1/1431720
for key2 in d2:
    d2[key2] *= 1/1342432
#   print("{:.3f}%".format(d1[key1]*100))
# print(d1)
# print(d2)
ds = [d1, d2]
d = {}
for k in d1.keys():
  d[k] = tuple(d[k] for d in ds)
# print(d)
print("Codon Frequency in Sal Frequency in Myc")
for i in d:
    print("{} {:.3f}% {:.3f}%".format(i,100*d[i][0],100*d[i][1]))
