#!/usr/bin/env python3

'''
Calcula estatisticas basicas dos tamanhos de sequencias de arquivos FASTA
Max, min, media e mediana

Entrada:    - Arquivo fasta (STDIN)
Saida:      - Estatisticas (STDOUT)
'''
import sys
from Bio import SeqIO
import statistics as st

fileName = sys.argv[1]
arraySize = []

for sequence in SeqIO.parse(fileName, 'fasta'):
    arraySize.append(len(sequence.seq))

print("Min\t{}".format(min(arraySize)))
print("Max\t{}".format(max(arraySize)))
print("Mean\t{}".format(st.mean(arraySize)))
print("Mode\t{}".format(st.mode(arraySize)))
print("Median\t{}".format(st.median(arraySize)))
print("Stdev\t{}".format(st.stdev(arraySize)))
#print(st.multimode(arraySize))
