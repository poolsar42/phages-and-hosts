#!/usr/bin/env python
# coding: utf-8

# In[17]:


import subprocess, os, argparse, pandas as pd, csv
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
import Bio
from Bio import UniGene

from Bio.Alphabet import generic_dna


# $\text{Searching bacteriophages in virus list and hosts Tax_ID list}$

# In[18]:


# Here I'm searching bacteriophages in the table (from Virus-Host DB FTP) and saving  result in list
# Also I saving host names in anpther list
z = []
i = 0
host_name = []
with open("virushostdb.tsv") as tsvfile:
    tsvreader = csv.reader(tsvfile, delimiter="\t",)
    for line in tsvreader:
        i += 1
        if ('Archaea' in line[9]):
            z.append(line)
            host_name.append(line[8])
        if ('Bacteria' in line[9]):
            z.append(line)
            host_name.append(line[8])
        if (line[9] == 'host lineage'):
            z.append(line)


# $\text{Searchinh names and genomes of bacheriophages}$

# $\text{First Method}$

# In[19]:


#Here I'm searching bacteriophages genomes using an FTP from Virus-Host DB
f = open('virushostdb.formatted.genomic.fna', 'r')
s = f.read()
f.close()
lst1 = []
lst2 = []
Genome = s.split('>')
Genome.remove('')
lst = []
for i in range(len(Genome)):
    sf = Genome[i].split('|')
    lst.append(sf)
for i in range(0,len(lst)):
    sf1 = lst[i][0].split(' ')
    lst1.append(sf1)
for i in range(0,len(lst1)):
    s = ''
    for j in range(1,len(lst1[i])):
        s += lst1[i][j] + ' '
    index = len(s)-1
    s = s[:index] + s[index+1:]
    lst2.append(s) #Virus names
l = 0
GENOMES = []
NAMES = []
for i in range(1,len(z)):
    p = lst2.index(z[i][1])
    NAMES.append(lst2[p]) #Bacteriophage names
    GENOMES.append(lst[p][5]) #Bacteriophage genomes


# $\text{Second Method (using biopython)}$

# In[20]:


Entrez.email = "sharafutdinov.ehn@phystech.edu"
Genes = []
#Can I do it in loop?
handle = Entrez.efetch(db="nucleotide", id=z[1][3], rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
Genes.append(repr(record.seq))
#print(Genes)


# $\text{Searching IDList from hosts names}$

# In[21]:


#I know host names and here I'm searching their genomes
Entrez.email = "sharafutdinov.ehn@phystech.edu"
HOST_IDs = []
s = ''
j = 0
#Can I do it in loop?
s = host_name[1] + "[Orgn]"
handle = Entrez.esearch(db="nucleotide", term=s, idtype="acc")
record = Entrez.read(handle)
record["Count"]
HOST_IDs.append(record["IdList"])
#print((HOST_IDs))


# $\text{I 'm looking for host genomes here}$

# In[22]:


Entrez.email = "sharafutdinov.ehn@phystech.edu"
Genes_Host = []
#for i in range(len(HOST_IDs)): And can I do this in loop?
handle = Entrez.efetch(db="nucleotide", id=HOST_IDs[0][1], rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()
Genes_Host.append(repr(record.seq))
#print((Genes_Host))

