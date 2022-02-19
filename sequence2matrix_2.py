#!/usr/bin/env python3

'''
Transforma as sequencias em um CSV de numeros
onde o ultimo numero eh o organismo, sendo
0 para mamiferos e 1 para bacteria e os demais
o valor da massa molecular de cada aminocido, 
a frequencia de cada aminoacido e os descritores
Z-class

Parametros: 
	Entrada: 	- arquivo multifasta
				- codir organismo (0 para mamiferos e 1 para bacterias)
	Saida:		- <STDOUT> aqruivo csv

Rodrigo Jardim
Maio 2021
'''
import sys
from Bio import SeqIO


# Aminoacidos com seus pesos moleculares e contagem de átomos (PUBCHEM)

#PubCHEM data
#Residue 	3-letter 	1-letter 	Molecular weight Heavy Atom Count
#https://pubchem.ncbi.nlm.nih.gov/compound/Glycine
#Glycine 	Gly 	G 	C2H3NO 	75.07 	5
#https://pubchem.ncbi.nlm.nih.gov/compound/Alanine
#Alanine 	Ala 	A 	C3H5NO 	89.09 	6
#https://pubchem.ncbi.nlm.nih.gov/compound/Serine
#Serine 	Ser 	S 	C3H5NO2 	105.09 	7
#https://pubchem.ncbi.nlm.nih.gov/compound/Proline
#Proline 	Pro 	P 	C5H7NO 	115.13 	8
#https://pubchem.ncbi.nlm.nih.gov/compound/Valine
#Valine 	Val 	V 	C5H9NO 	117.15 	8
#https://pubchem.ncbi.nlm.nih.gov/compound/l-Threonine
#Threonine 	Thr 	T 	C4H7NO2 	119.12 	8
#https://pubchem.ncbi.nlm.nih.gov/compound/Cysteine
#Cysteine 	Cys 	C 	C3H5NOS 	121.16 	7
#https://pubchem.ncbi.nlm.nih.gov/compound/l-Isoleucine
#Isoleucine 	Ile 	I 	C6H11NO 	131.17 	9
#https://pubchem.ncbi.nlm.nih.gov/compound/Leucine
#Leucine 	Leu 	L 	C6H11NO 	131.17 	9
#https://pubchem.ncbi.nlm.nih.gov/compound/Asparagine
#Asparagine 	Asn 	N 	C4H6N2O2 	132.12 	9
#https://pubchem.ncbi.nlm.nih.gov/compound/Aspartic-acid
#Aspartic acid 	Asp 	D 	C4H5NO3 	133.10 	9
#https://pubchem.ncbi.nlm.nih.gov/compound/Glutamine
#Glutamine 	Gln 	Q 	C5H8N2O2 	146.14 	10
#https://pubchem.ncbi.nlm.nih.gov/compound/Lysine
#Lysine 	Lys 	K 	C6H12N2O 	146.19 	10
#https://pubchem.ncbi.nlm.nih.gov/compound/Glutamic-acid
#Glutamic acid 	Glu 	E 	C5H7NO3 	147.13 	10
#https://pubchem.ncbi.nlm.nih.gov/compound/Methionine
#Methionine 	Met 	M 	C5H9NOS 	149.29 	9
#https://pubchem.ncbi.nlm.nih.gov/compound/Histidine
#Histidine 	His 	H 	C6H7N3O 	155.11 	11
#https://pubchem.ncbi.nlm.nih.gov/compound/Phenylalanine
#Phenylalanine 	Phe 	F 	C9H9NO 	165.19 	12
#https://pubchem.ncbi.nlm.nih.gov/compound/Arginine
#Arginine 	Arg 	R 	C6H12N4O 	174.20 	12
#https://pubchem.ncbi.nlm.nih.gov/compound/Tyrosine
#Tyrosine 	Tyr 	Y 	C9H9NO2 	181.19 	13
#https://pubchem.ncbi.nlm.nih.gov/compound/Tryptophan
#Tryptophan 	Trp 	W 	C11H10N2O 	204.22 	15 

# Molecular weight
dictAminoacids = {'g':75.07,'a':89.09, 's':105.09, 'p':115.13, 'v':117.15, 't':119.12, 'c':121.16,'i':131.17, 'l':131.17, 'n':132.12, 'd':133.10, 'q':146.14,'k': 146.19,'e':147.13, 'm':149.29, 'h':155.11, 'f':165.19, 'r':174.20, 'y':181.19, 'w':204.22, 'x':0, 'u':0, 'b':0}

#Heavy Atom Count
dictAminoacids_heavy_atom_count = {'g':5,'a':6, 's':7, 'p':8, 'v':8, 't':8, 'c':7,'i':9, 'l':9, 'n':9, 'd':9, 'q':10, 'k': 10, 'e':10, 'm':9, 'h':11, 'f':12, 'r':12, 'y':13, 'w':15, 'x':0, 'u':0, 'b':0}

# Aminoacidos com suas massas medias (Average mass)
# http://www2.riken.jp/BiomolChar/Aminoacidmolecularmasses.htm
#dictAminoacids = {'g':57.051,'a':71.078,'l':113.158,'m':131.196,'f':147.174,'w':186.210,'k':128.172,'q':128.129,'e':129.114,'s':87.077,'p':97.115,'v':99.131,'i':113.158,'c':103.143,'y':163.173,'h':137.139,'r':156.186,'n':114.103,'d':115.087,'t':101.104, 'x':0, 'u':0, 'b':0}

# Colunas da Zclass: ‘hydrophobicity’, ‘bulk of side chain’, ‘electronic properties’ 
# Artigo na pasta docs
dictZclass = {
	'f':[-4.92,1.3,0.45],
	'w':[-4.75,3.65,0.85],
	'i':[-4.44,-1.68,-1.03],
	'l':[-4.19,-1.03,-0.98],
	'v':[-2.69,-2.53,-1.29],
	'm':[-2.49,-0.27,-0.41],
	'y':[-1.39,2.32,0.01],
	'p':[-1.22,0.88,2.23],
	'a':[0.07,-1.73,0.09],
	'c':[0.71,-0.97,4.13],
	't':[0.92,-2.09,-1.4],
	's':[1.96,-1.63,0.57],
	'q':[2.19,0.53,-1.14],
	'g':[2.23,-5.36,0.3],
	'h':[2.41,1.74,1.11],
	'k':[2.84,1.41,-3.14],
	'r':[2.88,2.52,-3.44],
	'e':[3.08,0.039,-0.07],
	'n':[3.22,1.45,0.84],
	'd':[3.64,1.13,2.36]
}

aminoacids = ['g','a','l','m','f','w','k','q','e','s','p','v','i','c','y','h','r','n','d','t']

# Aminoacidos que nao existem mas que estao nas sequencias
#aminoacids.append('x')
#aminoacids.append('u')
#aminoacids.append('b')

import sys
import hashlib
from Bio import SeqIO

strFileIn = sys.argv[1]		# Arquivo multifasta em aminoacido
strOrganism = sys.argv[2]   # 0 para bacteria e 1 para mamifero

# Zerando a frequencia de cada aminoacido
frequence = {}
for letter in aminoacids:
	frequence[letter] = 0

# Montando a saida
output = {}
output['hash'] = []
output['descriptors'] = {}

# Cabecalho com as letras dos aminoacidos
strCSV = ''
for letter in aminoacids:
	strCSV+=letter
	strCSV+=","
	output['descriptors'][letter] = []

# Preenchendo a matriz
for sequence in SeqIO.parse(strFileIn, "fasta"):
	seq = str(sequence.seq).lower()
	
	# Calculando a frequencia de cada aminoacido
	for letter in seq:
		try:
			frequence[letter]+=1
		except:
			pass

	# Calculando o hash para a sequencia
	hash = hashlib.md5(seq.encode()).hexdigest()

	# Monta vetor de saida
	for letter in aminoacids:
		line = []
		
		# Inclui uma coluna com o valor ponderado de cada aminoacido
		# pelo seu peso molecular (dictAminoacids)
		line.append(str(dictAminoacids[letter] * frequence[letter]))
		
		# Inclui 3 colunas com os descritores da tabela Z, ponderado
		# pela frequencia de cada aminoacido
		for z in dictZclass[letter]:
			line.append(str(z * frequence[letter]))

		output['descriptors'][letter].append(line)
	
	output['hash'].append(hash)

# Imprimindo a saida
max = len(output['descriptors']['g'])
for i in range(max):
	for k, v in output.items():
		if (k=='hash'):
			# Nao estou imprimindo o hash da sequencia
			pass
			#print(output['hash'][i], end = ' ')
		else:
			for letter in aminoacids:
				lista = output['descriptors'][letter][i]
				print(','.join(output['descriptors'][letter][i]), end=',')

	print(strOrganism)