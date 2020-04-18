def read_sequences(file_path):
    sequence=[]
    file=open(file_path,'r')
    data=file.readlines()
    for line in data:
        if(line[0]!='>'):
            line=line[:-1]
            sequence.append(line)
    file.close()
    return sequence

import random

def background_frequency(pos,seq_dict,sequence,motif_length,background_count):
    for j in range(len(sequence)):
        if(j>=pos and j<=(pos+motif_length-1)):
            pass
        else:
            background_count=background_count+1
            temp_seq=sequence[j]+'0'
            try:
                seq_dict[temp_seq]=seq_dict[temp_seq]+1
            except:
                seq_dict[temp_seq]=1
    return background_count

def build_weight_matrix(selected_sequence,motif_length,sequence_motif_dict):
    seq_dict={}
    background_count=0
    for index in range(len(sequence)):
        if(index!=selected_sequence):
            if(len(sequence_motif_dict[index])==0):
                pos=random.randint(0,len(sequence[index])-motif_length)
            else:
                pos=sequence_motif_dict[index].index(max(sequence_motif_dict[index]))
            background_count+=background_frequency(pos,seq_dict,sequence[index],motif_length,0)
            for j in range(motif_length):
                temp_seq=sequence[index][j+pos]+str(j+1)
                try:
                    seq_dict[temp_seq]=seq_dict[temp_seq]+1
                except:
                    seq_dict[temp_seq]=1
    return background_count,seq_dict


def calculate_probability(n,b,B,seq_dict,background_count):
    g=['A','T','C','G']
    for index in range(7):
        for j in g:
            try:
                seq_dict[j+str(index+1)]=(seq_dict[j+str(index+1)]+b)/(n-1+B)
            except:
                seq_dict[j+str(index+1)]=(b/(n-1+B))
    for j in g:
        seq_dict[j+'0']=(seq_dict[j+'0']+b)/(background_count+B)  
    return seq_dict

def estimate_motif(prob_dict,sequence_index,motif_length):
    motif_list=[]
    for j in range(len(sequence[sequence_index])-motif_length):
        num=1
        den=1
        for index in range(motif_length):
            temp_seq=sequence[sequence_index][j+index]+str(index+1)
            try:
                num=num*prob_dict[temp_seq]
            except:
                num=num*(0.2/(76))
            temp_seq=sequence[sequence_index][j+index]+'0'
            den=den*prob_dict[temp_seq]
        motif_list.append(num/den)
    return motif_list

def normalize_list(motif_list):
    normalize=0
    for j in motif_list:
        normalize=normalize+j
    temp_motif_list=[]
    for j in motif_list:
        temp_motif_list.append(j/normalize)
    return temp_motif_list

import math
    

def find_motif_score(sequence_motif_dict):
    temp_score=0
    for motif in sequence_motif_dict:
        temp_score=temp_score+math.log2(max(sequence_motif_dict[motif]))
        
    print(temp_score)

def print_motif(prob_dict,motif_length):
    g=['A','T','C','G']
    motif=""
    for index in range(motif_length):
        max_score=0
        for j in g:
            temp_index=j+str(index)
            if(prob_dict[temp_index]>max_score):
                temp_motif=j
                max_score=prob_dict[temp_index]
        motif=motif+temp_motif
    print(motif)
            
         
def print_sequence_motifs(sequence,sequence_motif_dict):
    for j in range(len(sequence)):
        pos=sequence_motif_dict[j].index(max(sequence_motif_dict[j]))
        print("gene->",j)
        print(sequence[j][pos:(pos+7)])
    
    
    
    
def motif_find():
    sequence=read_sequences('C:\\Users\\ADMIN\\Downloads\\seqsA.fa')
    motif_dict={}
    sequence_motif_dict={}
    for j in range(len(sequence)):
        sequence_motif_dict[j]=[]
    for iteration in range(1000):
        j=random.randint(0,len(sequence)-1)
        background_count,seq_dict=build_weight_matrix(j,7,sequence_motif_dict)
        prob_dict=calculate_probability(75,0.5,2,seq_dict,background_count)
        motif_list=estimate_motif(prob_dict,j,7)
        sequence_motif_dict[j]=motif_list
        motif_list=normalize_list(motif_list)
        max_motif=max(motif_list)
        max_index=motif_list.index(max_motif)
        motif=str(j+1)+" "+str(max_index)
        try:
            motif_dict[motif]=motif_dict[motif]+1
        except:
            motif_dict[motif]=1
    print(prob_dict)
    find_motif_score(sequence_motif_dict)
    print_motif(prob_dict,7)
    print("")
    print("ALL sequence motifs")
    print_sequence_motifs(sequence,sequence_motif_dict)
    
            
motif_find()

