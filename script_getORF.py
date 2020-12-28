# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:02:47 2020

@author: Carolina Albuquerque Massena Ribeiro
"""


#!/usr/bin/env python3
import sys #importamos o modulo sys

file = sys.argv[1]

with open(file, 'r') as file:
    fasta = file.read() #Abrimos o arquivo .fasta

sequences = []
sequence_names = []

fasta = fasta.split('\n')

current_sequence = ''
for i, sequence in enumerate(fasta): #loop for
    if '>' in sequence and i>0:
        sequences.append(current_sequence)
        sequence_names.append(sequence)
        current_sequence = ''
    elif '>' not in sequence:
        current_sequence += sequence

end_seq = ['TAG', 'TAA', 'TGA'] #Sequencias de parada

real_sequences = []
index_pair = []
for sequence in sequences:
    for i in range(len(sequence)-2):
        current_start_seq = sequence[i:i+3]
        if current_start_seq == 'ATG': #Sequencia de inicio
            start_index = i
            end_index = None
            for j in range(i+3, len(sequence)-2, 3):     
                current_end_seq = sequence[j:j+3]
                if current_end_seq in end_seq:
                    end_index = j + 3
                    break                        
                        
            if end_index is not None:
                index_pair.append([start_index, end_index])
                undivided_sequence = sequence[start_index:end_index]
                divided_sequences = []
                
                for k in range(0, len(undivided_sequence)-2, 3):
                    divided_sequences.append(undivided_sequence[k:k+3])
                    
                real_sequences.append(divided_sequences)
                break        

#Fazendo um dicionario para parear as trincas com seu respectivo aminoacidos
trans_dict = {'AUG':'M', 'UAG, UAA, UGA':'-','UUU, UUC':'F', 'UUA, UUG, CUU, CUC, CUA, CUG':'L', 'AUU, AUC, AUA':'I', 'GUU, GUC, GUA, GUG':'V', 'UCU, UCC, UCA, UCG':'S', 'CCU, CCC, CCA, CCG':'P', 'ACU, ACC, ACA, ACG':'T', 'GCY, GCC, GCA, GCG':'A', 'UAU, UAC':'Y', 'CAU, CAC':'H', 'CAA, CAG':'Q', 'AAU, AAC':'N', 'AAA, AAG':'K', 'GAU, GAC':'D', 'GAA, GAG':'E', 'UGU, UGC':'C', 'UGG':'W', 'CGU, CGC, CGA, CGG, AGA, AGG':'R', 'AGU, AGC':'S', 'GGU, GGC, GGA, GGG':'G'}
big_len = 0
for sequence, index_pair, name in zip(real_sequences, index_pair, sequence_names): #loop for para determinar a maior das sequencias em nosso codigo
    if len(sequence) > big_len:
        big_len = len(sequence)
        longest_sequence = sequence
        se_index_pair = index_pair
        sequence_name = name
        

longest_sequence_trans = [] #loop for para substituir a base T, para a base U. Em nosso dicionario usamos as trincas com U
for trinca in longest_sequence:
    new_trinca = trinca.replace('T', 'U')
    for key, value in trans_dict.items():
        if new_trinca in key:
            longest_sequence_trans.append(value)
            

start = str(int(se_index_pair[0]) + 1)
end = se_index_pair[1]
frame_num = len(real_sequences)

#Criando os arquivos de texto para salvar as informações
        
with open('ORF.faa', 'a') as new_file: 
    new_file.write(f"{sequence_name}_frame{frame_num}-{start}-{end}\n")
    for trinca_trans in longest_sequence_trans:
        new_file.write(trinca_trans)

with open('ORF.fna', 'a') as new_file2:
    new_file2.write(f"{sequence_name}_frame{frame_num}-{start}-{end}\n")
    for trinca in longest_sequence:
        new_file2.write(trinca)

#Um print no console para informar o usuario que o codigo acabou       
print('done')