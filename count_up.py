#!/usr/bin/env python3

# download the gff file
import gzip, os, urllib.request
loadfile = "ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
gff = urllib.request.urlopen(loadfile)

# this is code which will count the number and length of genes in gff
total_gene_length = 0
with gzip.open(gff) as f:
        i=0
        for line in f:
                linestrip = line.decode('UTF-8').strip().split('\t')
                if linestrip[0][0] == "#":
                        continue
                elif linestrip[2] == "gene":
                        i += 1
                        genelen = abs(int(linestrip[4])-int(linestrip[3]))
                        total_gene_length += genelen
        print("total gene number in the gff is {}".format(i))
        print("Total length of genes is {}".format(total_gene_length))
# download the fasta file
fastafile="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"
fasta = urllib.request.urlopen(fastafile)
genome = ""
with gzip.open(fasta) as fa:
    for line in fa:
        if line[0] != ">":
            line = line.strip().upper()
            l = str(line)
            genome += l
print("The length of genome in fasta file is {}".format(len(genome)))
# compute the CDS length
total_CDS_length = 0
with gzip.open(gff) as f:
        i=0
        for line in f:
                linestrip = line.decode('UTF-8').strip().split('\t')
                if linestrip[0][0] == "#":
                        continue
                elif linestrip[2] == "CDS":
                        i += 1
                        CDSlen = abs(int(linestrip[4])-int(linestrip[3]))
                        total_CDS_length += CDSlen
# cumpute the percentage of genome which is coding
percentage = 100 * total_CDS_length/len(genome)
p = round(percentage, 3)
print("The percentage of genome which is coding is {}%".format(p, '.3f'))
