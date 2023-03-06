# Introduction
This is Homework 1 for S&DS 352: Biomedical Data Mining at Yale University, Spring 2023.

Written by Danielle Lee (MY '24)

This is a naive implementation of the Smith-Waterman algorithm for local alignment of two protein sequences. It supports an affine gap penalty and a similarity matrix.

The main script takes in two required and two optional arguments, for a total of four maximum arguments. The required inputs are the input file containing the two protein sequences, and the score file containing the similarity matrix for all amino acid pairings. The optional inputs are the gap opening penalty and gap extension penalty, which should be negative integers. When these inputs are not provided, the defaults are gap opening penalty = -2 and gap extension penalty = -1.

# Usage
```
$ python smith_waterman.py -i <input file> -s <score file> [-o <open gap penalty>] [-e <ext gap penalty]
```

# Example
Using the input and score files in this repo,
```
python smith_waterman.py -i input.txt -s blosum62.txt
```
Running the above line should yield a file called out.txt, the last few lines of which should look like this:
```
----------------------
|Best Local Alignment|
----------------------
Alignment Score:135
Alignment Results:
    (SSSP--QP---KK-KPL---DGEYFTLQIRGRERFEMFRELNEALELKD---AQAGKEPGGSRAHS-S)
     || |  ||   || |       | |||| ||  ||||    |  |||||   | |         || |
TFKE(SSLPSSQPEDSKKTKTTSSNEEEIFTLQVRGKKRFEMLKMINDSLELKDLVPA-ADQDKYRQKLHSKS)TS
```

# Installation
```
git clone https://github.com/leedanielle/smith_waterman
```

