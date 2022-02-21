#!/usr/bin/env python3

'''
Filtra as sequencias de arquivo FASTA por tamanho das sequencias

Entrada:    - Nome do arquivo fasta 
            - media
            - desvio padrao em percentual

Saida:      - arquivo fasta filtrado pelos tamanhos (menor e maior) (STDOUT)
'''
import sys
from Bio import SeqIO

fileName = sys.argv[1]
mean = float(sys.argv[2])
stdev = float(sys.argv[3])

minSize = mean - (mean * stdev / 100)
maxSize = mean + (mean * stdev / 100)
for sequence in SeqIO.parse(fileName, 'fasta'):
    if (len(sequence.seq) >= minSize) and (len(sequence.seq) <= maxSize):
        print(">{}{}\n{}".format(sequence.id, sequence.description, sequence.seq))
