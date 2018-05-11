# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 19:12:53 2017

@author: TatlierY
"""


#%%
# Functions for Project

def build_scoring_matrix(alphabet, diag_score, off_diag, dash_score):
    """
    Takes as input a set of characters alphabet and three scores diag_score, off_diag_score, and dash_score. 
    The function returns a dictionary of dictionaries whose entries are indexed by pairs of characters in 
    alphabet plus '-'. The score for any entry indexed by one or more dashes is dash_score. 
    The score for the remaining diagonal entries is diag_score. Finally, the score for the remaining off-diagonal 
    entries is off_diag_score.
    """
    
    scoring_dict={}
    dash='-'
    alphabet_copy=set(alphabet)
    alphabet_copy.update(dash)
    score=None    
    
    for character_i in alphabet_copy:
        scoring_dict[character_i]={}
        for character_j in alphabet_copy:
            if dash in [character_i,character_j]:
                score=dash_score
            else:
                if character_i==character_j:
                    score=diag_score
                else:
                    score=off_diag

            scoring_dict[character_i][character_j]=score
    
    return scoring_dict        

def check_flag(value,global_flag):
    """
    Makes a distinction between two different methods (Q8 and Q12)
    """
    if(global_flag==True):
        return value
    elif(global_flag==False):
        return max(0,value)
    else:
        print("global flag value is not allowed")
            
            
def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag): 
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix 
    scoring_matrix. The function computes and returns the alignment matrix for seq_x and seq_y as described in 
    the Homework. If global_flag is True, each entry of the alignment matrix is computed using the method 
    described in Question 8 of the Homework. If global_flag is False, each entry is computed using the method 
    described in Question 12 of the Homework.
    """
    row_len, col_len = len(seq_x), len(seq_y)
    alignment_matrix=[[0 for _ in range(col_len+1)] for _ in range(row_len+1)]
    dash='-'
    
    #Pass the elements for the first column
    for index_i in range(1,row_len+1):
        alignment_matrix[index_i][0]=check_flag(alignment_matrix[index_i-1][0]+scoring_matrix[seq_x[index_i-1]][dash],global_flag) 
    
    #Pass the elements for the first row
    for index_j in range(1,col_len+1):
        alignment_matrix[0][index_j]=check_flag(alignment_matrix[0][index_j-1]+scoring_matrix[dash][seq_y[index_j-1]],global_flag) 
    
    #Filling all other elements    
    for index_i in range(1,row_len+1):
        for index_j in range(1,col_len+1):
            element_ij=check_flag(alignment_matrix[index_i-1][index_j-1]+scoring_matrix[seq_x[index_i-1]][seq_y[index_j-1]],global_flag)
            element_ij_10=check_flag(alignment_matrix[index_i-1][index_j]+scoring_matrix[seq_x[index_i-1]][dash],global_flag)
            element_ij_01=check_flag(alignment_matrix[index_i][index_j-1]+scoring_matrix[dash][seq_y[index_j-1]],global_flag)
            alignment_matrix[index_i][index_j]=max(element_ij,element_ij_10,element_ij_01)
            
    return alignment_matrix

    
def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix): 
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix. This function computes 
    a global alignment of seq_x and seq_y using the global alignment matrix alignment_matrix.The function returns a tuple of the form (score, align_x, align_y) 
    where score is the score of the global alignment align_x and align_y. Note that align_x and align_y should have the same length and may include the 
    padding character '-'.
    """
    row_len, col_len = len(seq_x), len(seq_y)
    x_ap, y_ap, dash ='','','-'  
    index_i, index_j = int(row_len), int(col_len)
    
    while index_i>0 and index_j>0:
        if(alignment_matrix[index_i][index_j]==alignment_matrix[index_i-1][index_j-1]+scoring_matrix[seq_x[index_i-1]][seq_y[index_j-1]]):
            x_ap = seq_x[index_i-1] + x_ap
            y_ap = seq_y[index_j-1] + y_ap
            index_i, index_j = index_i-1, index_j-1
        else:
            if(alignment_matrix[index_i][index_j]==alignment_matrix[index_i-1][index_j]+scoring_matrix[seq_x[index_i-1]][dash]):
                x_ap=seq_x[index_i-1]+x_ap
                y_ap=dash+y_ap
                index_i -= 1
            else:
                x_ap=dash+x_ap
                y_ap = seq_y[index_j-1]+y_ap
                index_j -= 1
    
    while index_i>0:
        x_ap, y_ap=seq_x[index_i-1]+x_ap, dash+y_ap
        index_i -= 1
        
    while index_j>0:
        x_ap, y_ap= dash+x_ap, seq_y[index_j-1]+y_ap
        index_j -= 1
    
    return (alignment_matrix[row_len][col_len],x_ap,y_ap)


#%%

def find_opt_cell(row_len, col_len, alignment_matrix):
    """
    Standard search algorithm in grid
    """
    opt_value, opt_cell=0, (0,0)
    
    for index_i in range(row_len+1):
        for index_j in range(col_len+1):
            if(alignment_matrix[index_i][index_j]>opt_value):
                opt_cell=(index_i,index_j)
                opt_value=alignment_matrix[index_i][index_j]
    
    return opt_value,opt_cell


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix): 
    """
    Takes as input two sequences seq_x and seq_y whose elements share a common alphabet with the scoring matrix scoring_matrix. 
    This function computes a local alignment of seq_x and seq_y using the local alignment matrix alignment_matrix.The function returns a tuple 
    of the form (score, align_x, align_y) where score is the score of the optimal local alignment align_x and align_y. Note that align_x and 
    align_y should have the same length and may include the padding character '-'.
    """
    
    row_len, col_len = len(seq_x), len(seq_y)
    x_ap, y_ap, dash ='','','-' 
    
    opt_value,opt_cell=find_opt_cell(row_len, col_len, alignment_matrix)
    index_i, index_j = opt_cell[0], opt_cell[1]
    
    while alignment_matrix[index_i][index_j] != 0:
        if(alignment_matrix[index_i][index_j]==alignment_matrix[index_i-1][index_j-1]+scoring_matrix[seq_x[index_i-1]][seq_y[index_j-1]]):
            x_ap = seq_x[index_i-1] + x_ap
            y_ap = seq_y[index_j-1] + y_ap
            index_i, index_j = index_i-1, index_j-1
        elif(alignment_matrix[index_i][index_j]==alignment_matrix[index_i-1][index_j]+scoring_matrix[seq_x[index_i-1]][dash]):
            x_ap=seq_x[index_i-1]+x_ap
            y_ap=dash+y_ap
            index_i -= 1
        else:
            x_ap=dash+x_ap
            y_ap=seq_y[index_j-1]+y_ap
            index_j -= 1
    
    return opt_value,x_ap,y_ap

#%%

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urllib.request.urlopen(filename)
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        line=str(line)[2:-5]
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][str(ykey)[2:-1]] = int(val)
    return scoring_dict

#Q1        
import urllib
from load_protein import PAM50_URL, HUMAN_EYELESS_URL, FRUITFLY_EYELESS_URL, CONSENSUS_PAX_URL,WORD_LIST_URL

alphabet={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','-'} 
#Read scoring matrix
scoring_matrix=read_scoring_matrix(PAM50_URL)

# load HUMAN_EYELESS_URL
HEU=urllib.request.urlopen(HUMAN_EYELESS_URL)
HEU_seq = str(HEU.read())
HEU_seq=HEU_seq[2:-5]

# load HUMAN_EYELESS_URL
FEU=urllib.request.urlopen(FRUITFLY_EYELESS_URL)
FEU_seq = str(FEU.read())
FEU_seq=FEU_seq[2:-5]

#
alignment_matrix_HEU_FEU=compute_alignment_matrix(HEU_seq, FEU_seq, scoring_matrix, False)
(score_HEU_FEU,cla1_HEU_FEU,cla2_HEU_FEU)=compute_local_alignment(HEU_seq, FEU_seq, scoring_matrix, alignment_matrix_HEU_FEU)

#Q2

# Load CONSENSUS_PAX_URL
CPU=urllib.request.urlopen(CONSENSUS_PAX_URL)
CPU_seq = str(CPU.read())
CPU_seq=CPU_seq[2:-1]

#Find dashes
dash='-'
dash_indices_cla1=[i for i,j in enumerate(cla1_HEU_FEU) if j==dash]
dash_indices_cla2=[i for i,j in enumerate(cla2_HEU_FEU) if j==dash]

#Better Approach, use replace method
cla1_cleaned=cla1_HEU_FEU.replace(dash, "")
cla2_cleaned=cla2_HEU_FEU.replace(dash, "")

#
alignment_matrix_CPU1=compute_alignment_matrix(CPU_seq,cla1_cleaned,scoring_matrix, True)
am_CPU_1,x_ap1,y_ap1=compute_global_alignment(CPU_seq, cla1_cleaned, scoring_matrix, alignment_matrix_CPU1)

#
alignment_matrix_CPU2=compute_alignment_matrix(CPU_seq,cla2_cleaned,scoring_matrix, True)
am_CPU_2,x_ap2,y_ap2=compute_global_alignment(CPU_seq, cla2_cleaned, scoring_matrix, alignment_matrix_CPU2)

def similarity(seq1,seq2):
    """
    computes similarity between two string (with same size)
    """
    diff_count=0
    
    for i in range(len(seq1)):
        if(seq1[i] != seq2[i]):
            diff_count += 1
    
    return 1-(diff_count/len(seq1))
        
#%%

#Q4-Q6
import random
import matplotlib.pyplot as plt

def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials,normalize=False):
    """
    Function takes as input two sequences seq_x and seq_y, a scoring matrix scoring_matrix, and a number of trials num_trials. 
    This function returns a dictionary scoring_distribution that represents an un-normalized distribution generated by performing the following process 
    num_trials times:
    1. Generate a random permutation rand_y of the sequence seq_y using random.shuffle().
    2. Compute the maximum value score for the local alignment of seq_x and rand_y using the score matrix scoring_matrix.
    3. Increment the entry score in the dictionary scoring_distribution by one.
    """
    trial_dict={}
    if(normalize==False):
        norm_factor=1
    else:
        norm_factor=num_trials
    
    for index_i in range(num_trials):
        print(index_i)
        perm_list=list(seq_y)
        random.shuffle(perm_list)
        perm_seq="".join(perm_list)
        alignment_matrix=compute_alignment_matrix(seq_x, perm_seq, scoring_matrix, False)
        opt_score,x_ap,y_ap=compute_local_alignment(seq_x, perm_seq, scoring_matrix, alignment_matrix)
        if(opt_score in trial_dict):
            trial_dict[opt_score] += 1/norm_factor
        else:
            trial_dict[opt_score] = 1/norm_factor
    
    return trial_dict
    
num_trials=1000
generate_dict_unnormalized=generate_null_distribution(HEU_seq,FEU_seq,scoring_matrix,num_trials)
generate_dict_normalized={}
for x,y in generate_dict_unnormalized.items():
    generate_dict_normalized[x]=y/num_trials

#generate_dict_normalized=generate_null_distribution(HEU_seq,FEU_seq,scoring_matrix,num_trials,True)
mean=sum([x*y/num_trials for x,y in generate_dict_unnormalized.items()])
var=sum([(y*(x-mean)**2)/num_trials for x,y in generate_dict_unnormalized.items()])
std=var**(1/2)
z_score=(875-mean)/std

#Make plot
hist_list=[]
for key,value in generate_dict_unnormalized.items():
    for index_i in range(int(value)):
        hist_list.append(key)
        
plt.hist(hist_list,bins=25,normed=True)

plt.xlabel('Score')
plt.ylabel('Normalized frequency')
plt.title('Normalized distribution cla for permutated strings (n=1000)')
plt.show()

#%%
#Q8

#Load data
WLU=urllib.request.urlopen("http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt")
WLU_read=WLU.read()
WLU_list=str(WLU_read).split('\\n')
WLU_list[0]="aa"
del WLU_list[-1]

def check_spelling_one(ref_word,word_list):
    """
    iterates through word_list and returns the set of all words that are within one edit distance dist of the string checked_word.
    """
    candidate_list=[]
    one_distance_list=[]
    alphabet=[]
    ref_len=len(ref_word)
    comparison_word=" "+str(ref_word)+" "
    
    #Construct a smaller search list for efficiency 
    for word in word_list:
        if(len(word) in [ref_len,ref_len+1,ref_len+2]):
            candidate_list.append(word)
    
    #Construct an alphabet for searches
    alphabet = []
    for letter in range(97,123):
        alphabet.append(chr(letter))
    
    for index_i in range(ref_len+2):
        for letter in [x for x in alphabet if x !=comparison_word[index_i]]:
            if(index_i==0):
                new_word=comparison_word[:index_i]+letter+comparison_word[index_i+1:]
            else:
                new_word=ref_word[:index_i-1]+letter+ref_word[index_i:]
                
            if(new_word in candidate_list):
                one_distance_list.append(new_word)

    return one_distance_list

def check_spelling_two(ref_word,word_list):
    """
    iterates through word_list and returns the set of all words that are within two edit distance dist of the string checked_word.
    """
    candidate_list=[]
    one_distance_list=[]
    alphabet=[]
    ref_len=len(ref_word)
    comparison_word=" "+str(ref_word)+" "
    
    #Construct a smaller search list for efficiency 
    for word in word_list:
        if(len(word) in [ref_len,ref_len+1,ref_len+2]):
            candidate_list.append(word)
    
    #Construct an alphabet for searches
    alphabet = [""]
    for letter in range(97,123):
        alphabet.append(chr(letter))
    
    for index_i in range(ref_len+2):
        for letter in [x for x in alphabet if x !=comparison_word[index_i]]:
            if(index_i==0):
                new_word=comparison_word[:index_i]+letter+comparison_word[index_i+1:]
            else:
                new_word=ref_word[:index_i-1]+letter+ref_word[index_i:]
                
            if(new_word in candidate_list):
                one_distance_list.append(new_word)

    return one_distance_list