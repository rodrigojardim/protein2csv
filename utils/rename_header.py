#!/usr/bin/env python3

'''
Renomeia as sequencias de arquivo FASTA

Entrada:    - Nome do arquivo fasta 
            - termo 

Saida:      - arquivo fasta renomeado (STDOUT)
'''
import sys
from Bio import SeqIO

fileName = sys.argv[1]
termo = sys.argv[2]
i = 1

for sequence in SeqIO.parse(fileName, 'fasta'):
    print(">{} {}\n{}".format(termo + "_" + str(i), sequence.description, sequence.seq))
    i += 1
