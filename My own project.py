#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess, os, argparse, pandas as pd, csv
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez
import Bio
from Bio import UniGene

from Bio.Alphabet import generic_dna


# $\text{Searching bacteriophages in virus list and hosts Tax_ID list}$

# In[34]:


# Here I'm searching bacteriophages in the table and saving  result in list
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
print((z[157][3]))


# $\text{Searchinh names and genomes of bacheriophages}$

# $\text{First Method}$

# In[35]:


#Here I'm searching bacteriophages genomes
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
    lst2.append(s) #Названия вирусов под всеми индексами
l = 0
GENOMES = []
NAMES = []
for i in range(1,len(z)):
    p = lst2.index(z[i][1])
    NAMES.append(lst2[p]) #Имена бактериофагов
    GENOMES.append(lst[p][5]) #Геномы бактериофагов


# $\text{Second Method}$

# In[36]:


Entrez.email = "sharafutdinov.ehn@phystech.edu"
Genes = []
for i in range(1, 2):
    handle = Entrez.efetch(db="nucleotide", id=z[i][3], rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    Genes.append(repr(record.seq))
print((Genes))


# $ \text{Seq test} $

# In[37]:


print(Bio.__version__)
my_seq = Seq("AGTACACTGGT")
my_seq


# $\text{Searching IDList from hosts names}$

# In[44]:


Entrez.email = "sharafutdinov.ehn@phystech.edu"
HOST_IDs = []
HOST_ID = []
s = ''
j = 0
for i in range(1,2):
    s = host_name[i] + "[Orgn]"
    handle = Entrez.esearch(db="nucleotide", term=s, idtype="acc")
    record = Entrez.read(handle)
    record["Count"]
    HOST_IDs.append(record["IdList"])
    HOST_ID.append(HOST_IDs[j][1])
    j += 1
print((HOST_ID))


# In[45]:


Entrez.email = "sharafutdinov.ehn@phystech.edu"
Genes_Host = []
for i in range(len(HOST_ID)):
    handle = Entrez.efetch(db="nucleotide", id=HOST_ID[i], rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    Genes_Host.append(repr(record.seq))
print((Genes_Host))


# $\text{Genes Searching}$
