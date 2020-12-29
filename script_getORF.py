# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 16:02:47 2020

@author: Carolina Albuquerque Massena Ribeiro
"""


#!/usr/bin/env python3
import sys # importamos o modulo sys

file = 'teste.fasta' # selecionando arquivo pela command line

with open(file, 'r') as file:
    fasta = file.read() # lendo o arquivo selecionado

sequences = [] # criando uma lista vazia aonde serão colocadas as sequencias do arquivo fasta
sequence_names = [] # criando uma lista vazia aonde serão coloados os nomes das sequencias do arquivo fasta

fasta = fasta.split('\n') #Divindo o string extraido do arquivo

#O bloco de código abaixo extrai os nomes das sequencias e as sequencias do arquivo fasta e popula as listas criadas anteriormente
current_sequence = ''
for i, sequence in enumerate(fasta): 
    if '>' in sequence and i>0:
        sequences.append(current_sequence)
        sequence_names.append(sequence)
        current_sequence = ''
    elif '>' not in sequence:
        current_sequence += sequence

end_seq = ['TAG', 'TAA', 'TGA'] #Códons de parada

# O bloco abaixo procura os trechos das sequencias que começam com 'ATG' e terminam com 'TAG' ou 'TAA' ou 'TGA'
real_sequences = [] # essa lista será populada com os trechos mencionados anteriormente
index_pair = [] # essa lista será populada com os indices de inicio e fim desses trechos
for sequence in sequences:
    i = 0
    current_valid_sequences = [] # aqui todas os trechos validos de uma sequencia serão colocados
    current_valid_index_pairs = [] # aqui todos os indices desses trechos serão colocados
    while i <= len(sequence)-3:
        current_start_seq = sequence[i:i+3]
        if current_start_seq == 'ATG':
            start_index = i
            end_index = None
            for j in range(i+3, len(sequence)-2, 3):     
                current_end_seq = sequence[j:j+3]
                if current_end_seq in end_seq:
                    end_index = j + 3
                    break                        
                        
            if end_index is not None:
                current_valid_index_pairs.append([start_index, end_index]) 
                undivided_sequence = sequence[start_index:end_index]
                divided_sequences = [] # aqui os trechos serão colocados dividos em trincas
                
                for k in range(0, len(undivided_sequence)-2, 3):
                    divided_sequences.append(undivided_sequence[k:k+3])
                    
                current_valid_sequences.append(divided_sequences)
                i = end_index
            else:
                i += 3
        else:
            i += 1
                
    real_sequences.append(current_valid_sequences)
    index_pair.append(current_valid_index_pairs)
    
#Fazendo um dicionario para parear as trincas com seu respectivo aminoacidos
trans_dict = {'AUG':'M', 'UAG, UAA, UGA':'-','UUU, UUC':'F', 'UUA, UUG, CUU, CUC, CUA, CUG':'L', 'AUU, AUC, AUA':'I', 'GUU, GUC, GUA, GUG':'V', 'UCU, UCC, UCA, UCG':'S', 'CCU, CCC, CCA, CCG':'P', 'ACU, ACC, ACA, ACG':'T', 'GCY, GCC, GCA, GCG':'A', 'UAU, UAC':'Y', 'CAU, CAC':'H', 'CAA, CAG':'Q', 'AAU, AAC':'N', 'AAA, AAG':'K', 'GAU, GAC':'D', 'GAA, GAG':'E', 'UGU, UGC':'C', 'UGG':'W', 'CGU, CGC, CGA, CGG, AGA, AGG':'R', 'AGU, AGC':'S', 'GGU, GGC, GGA, GGG':'G'}
# abaixo são criadas as listas onde serão guardadas as informações referentes ao trechos mais longos de cada sequencia
longest_stretches = []
longest_index_pairs = []
longest_sequence_names = []

# o código abaixo descobre qual o trecho mais longo de cada sequencia e armazena suas informações
for sequence, indexes, name in zip(real_sequences, index_pair, sequence_names):
    big_len = 0
    for stretch, index_pair in zip(sequence, indexes):
        if len(stretch) > big_len:
            big_len = len(stretch)
            longest_stretch = stretch
            se_index_pair = index_pair
            
    longest_stretches.append(longest_stretch) 
    longest_index_pairs.append(se_index_pair)
    longest_sequence_names.append(name)

# o código abaixo converte os trechos para os respectivos aminoacidos
longest_stretches_trans = []
for stretch in longest_stretches:
    stretch_trans = '' #loop for para substituir a base T, para a base U. Em nosso dicionario usamos as trincas com U
    for trinca in stretch:
        new_trinca = trinca.replace('T', 'U')
        for key, value in trans_dict.items():
            if new_trinca in key:
                stretch_trans += value
                
    longest_stretches_trans.append(stretch_trans)        

#Criando os arquivos de texto para salvar as informações
        
with open('ORF.faa', 'a') as new_file:
    for stretch_trans, sequence_name, index_pair in zip(longest_stretches_trans, longest_sequence_names, longest_index_pairs):
        start = index_pair[0]+1
        end = index_pair[1]
        frame_num = len(stretch_trans)
        new_file.write(f"{sequence_name}_frame{frame_num}-{start}-{end}\n")
        new_file.write(stretch_trans + '\n')

with open('ORF.fna', 'a') as new_file2:
    for stretch, sequence_name, index_pair in zip(longest_stretches, longest_sequence_names, longest_index_pairs):
    
        start = index_pair[0]+1
        end = index_pair[1]
        frame_num = len(stretch)
        
        new_file2.write(f"{sequence_name}_frame{frame_num}-{start}-{end}\n")
        for trinca in stretch:
            new_file2.write(trinca)
        new_file2.write('\n')
#Um print no console para informar o usuario que o codigo acabou       
print('done')