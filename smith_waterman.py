#!/usr/bin/python
__author__ = "Danielle Lee"
__email__ = "danielle.lee@yale.edu"
__copyright__ = "Copyright March 2023"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python DL_hw1.py -i <input file> -s <score file>
### Example: python DL_hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
import numpy as np
import itertools
from pandas import *
from contextlib import redirect_stdout

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

# Calculate the affine gap penalty
def gap_penalty(distance, opengap, extgap): 
    return opengap+(extgap*(distance-1))

# Calculate the score matrix and traceback.
# Optimal alignment is computed using traceback and printed to stdout.
def matrix(a, b, similarity, opengap, extgap):
    
    # a and b are the protein sequences.
    # Throughout this function we assume a is shorter than b.
    # Thus, switch a and b if a is longer than b.
    if len(a) > len(b):
        temp = a
        a = b
        b = temp
    # We assume the amino acids are universal.    
    amino_acids = 'ABCDEFGHIKLMNPQRSTVWXYZ'
    
    # Initialize score and pointer matrix
    score_matrix = np.zeros((len(a) + 1, len(b) + 1), int)
    
    # This will later tell us where to go during traceback
    pointer_matrix = np.zeros((len(a) + 1, len(b) + 1), int)  
    
    max_score = 0
    
    # Iterate through all amino acid pairings
    for i, j in itertools.product(range(1, score_matrix.shape[0]), range(1, score_matrix.shape[1])):
        
        # Find the match score between every pair and add to diagnoal
        a_i = amino_acids.index(a[i-1])
        b_j = amino_acids.index(b[j-1])
        match_score = similarity[a_i][b_j]
        match = score_matrix[i-1, j-1] + match_score
        
        # List of all match scores to the left accounting for gap penalty
        scores_left = [(score_matrix[i,j-l]+gap_penalty(l, opengap, extgap)) for l in range(1,i+1)]
        
        # Same but for all upwards scores
        scores_up = [(score_matrix[i-k,j]+gap_penalty(k, opengap, extgap)) for k in range(1,i+1)]
        
        # Find the leftwards and upwards maxima
        leftmax = max(scores_left)
        upmax = max(scores_up)
        
        score_matrix[i,j] = max(match, leftmax, upmax, 0)
        
        if score_matrix[i,j] == 0:
            pointer_matrix[i,j] = 0 # 0 means end of the trace
            
        if score_matrix[i,j] == upmax:
            pointer_matrix[i,j] = 1 # 1 means the score came from up
            
        if score_matrix[i,j] == leftmax:
            pointer_matrix[i,j] = 2 # 2 means the score came from left
            
        if score_matrix[i,j] == match:
            pointer_matrix[i,j] = 3 # 3 means score came from diagonal
        
        # Update max score
        if score_matrix[i,j] >= max_score:
            max_score = score_matrix[i,j]
            max_i, max_j = i, j
    
    # Now time to align the two using traceback
    alignA, alignB = '', ''
    i, j = max_i, max_j
    
    while pointer_matrix[i,j] != 0:
        if pointer_matrix[i,j] == 3: # Trace diagonally
            alignA += a[i-1]
            alignB += b[j-1]
            i -= 1
            j -= 1
        elif pointer_matrix[i,j] == 2: # Trace leftwards
            alignA += '-'
            alignB += b[j-1]
            j -= 1
        elif pointer_matrix[i,j] == 1: # Trace upwards
            alignA += a[i-1]
            alignB += '-'
            i -= 1
    
    # Reverse both alignments, since we've added the characters backwards
    alignA = alignA[::-1]
    alignB = alignB[::-1]

    gap = j-i
    
    # Now draw the lines between the two alignments
    align_output = ''
    alignment = ''
    for index in range(0, len(alignA)):
        if alignA[index] == alignB[index]:
            alignment += '|'
        elif alignA[index] != alignB[index]:
            alignment += ' '

    if gap > 0:
        align_output += ' '*np.absolute(gap)
        
    align_output += a[:i] + '(' + alignA + ')' + a[max_i:]
    align_output += '\n'
    align_output += ' '*(max(i,j)+1) + alignment
    align_output += '\n'
    
    if gap < 0:
        align_output += ' '*np.absolute(gap)

    align_output += b[:j] + '(' + alignB + ')' + b[max_j:]
    align_output += '\n'
    
    construct_output(a, b, score_matrix)
    print("Alignment Score:", end='')
    print(max_score)
    print("Alignment Results:")
    print(align_output)
    return

# Construct the texts of the output and the tab-delimited score matrix
def construct_output(a, b, score_matrix):
    print('-----------')
    print('|Sequences|')
    print('-----------')
    print('sequence1')
    print(a)
    print('sequence2')
    print(b)
    print('--------------')
    print('|Score Matrix|')
    print('--------------')

    for i in range(0, len(a)+1):
        if i-1 >= 0:
            print('\t' + a[i-1], end='')
    print()
    for j in range(0, len(b)+1):
        if j-1 >= 0:
            print(b[j-1], end='')
        for i in range(0, len(a)+1):
            print('\t', end='')
            print(score_matrix[i, j], end='')
        print()
    print('----------------------')
    print('|Best Local Alignment|')
    print('----------------------')

### Put everything together
def runSW(inputFile, scoreFile, openGap, extGap):
    
    # Open the input file and get the two sequences, seqA and seqB
    myfile = open(inputFile, 'r')
    lines = myfile.readlines()
    seqA = lines[0].strip()
    seqB = lines[1].strip()
    
    # Open the score file and read in the similarity matrix
    similarity = DataFrame(np.loadtxt(scoreFile,
        dtype='int', skiprows=1, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)))
    
    # Open outgoing file and write in the matrix
    with open('sw_aligned.txt', 'w') as f:
        with redirect_stdout(f):
            matrix(seqA, seqB, similarity, openGap, extGap)
    f.close()
    
    
### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)
