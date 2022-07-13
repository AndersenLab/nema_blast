##### Import Packages

# Standard math, numpy, sys, os, pandas
import math
import glob
import numpy as np
import sys
import os
import time
from datetime import timedelta, date
import pandas as pd

# From Biopython
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.substitution_matrices import Array


##### Global Variables
alignment_list = []
df = pd.DataFrame(columns = ["sequence", "best_match", 
    "sequence_length", "best_match_length", "score"])
divide_bool = False
check_reverse_complement = True
specific = False

# Set up the directory path
path = os.getcwd()
os.chdir(path)

# This section, specifically the reference_input_path variable, can be changed
# to adjust the default reference file or 
# ______________________________________________________________________________

# Pre-setup file paths
# query_input_path = ""
# reference_input_path = "reference.fa"

if len(sys.argv) > 1:
    # Get files from command line
    query_input_path = sys.argv[1]
    # reference_input_path = sys.argv[2]

    # Try-catch mechanism for determining reference_input_path
    try:
        secondarg = sys.argv[2]
        if secondarg != "-m":
            reference_input_path = sys.argv[2]
        else:
            reference_input_path = "bin/reference.fa"
    except:
        reference_input_path = "bin/reference.fa"
    
    try:
        modify = sys.argv[3]
        if modify == "-m":
            specific = True
        else:
            specific = False
    except:
        try:
            modify = sys.argv[2]
            if modify == "-m":
                specific = True
            else:
                specific = False
        except:
            specific = False

else:
    query_input_path = str(input("\nEnter path to query file/directory: "))
    reference_input_path = str(input("\nEnter path to query file/directory: "))
    modify = str(input("\nManually modify scoring system? [y/n]: "))
    if modify == "y":
        specific = True
    else:
        specific = False

# ______________________________________________________________________________

print("\n\nQuery path received: " + str(os.path.basename(query_input_path)))
print("Reference path received: " + str(os.path.basename(reference_input_path)))


##### Initialize the PairwiseAligner object
aligner = Align.PairwiseAligner()
dnafull_nums_array = np.array([
        5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, 
        -4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 
        -4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2, 
        -4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, 
        -4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, 
        1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, 
        1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, 
        -4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, 
        -4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, 
        1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, 
        -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, 
        -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, 
        -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, 
        -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, 
        -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
        ])
dnafull = Array("ATGCSWRYKMBVHDN", 2, dnafull_nums_array.reshape(15, 15))
'''
dnafull (also known as EDNAFULL) is a scoring matrix used by EMBOSS needle
alignment tool
'''

def setUpAligner():
    '''
    This function sets up the PairwiseAligner object based on user input.

    There are two default settings: 
    1. EMBOSS Needle scoring system from 
        https://www.ebi.ac.uk/Tools/psa/emboss_needle/
    2. blastn default scoring system as detailed in
        https://biopython.org/DIST/docs/tutorial/Tutorial.html

    See the Biopython documentation above to understand gap scoring. 
    Accurate alignment is dependent on choosing the correct scoring system.
    '''
    global check_reverse_complement
    global divide_bool
    global specific
    
    if specific:
        revcomp = str(input("\nConsider reverse complements? [y/n]: "))
        if revcomp == "y":
            check_reverse_complement = True
        else:
            check_reverse_complement = False
        # check_reverse_complement is global and will be used later with the 
        # scoring system

        needle = str(input("\nUse EMBOSS Needle default scoring system? [y/n]: "))
        if needle == "y":
            aligner.substitution_matrix = dnafull
            aligner.internal_open_gap_score = -10
            aligner.internal_extend_gap_score = -0.5
            aligner.left_open_gap_score = 0
            aligner.left_extend_gap_score = 0
            aligner.right_open_gap_score = 0
            aligner.right_extend_gap_score = 0
            print("\nAlignment algorithm:", aligner.algorithm)
            return

        blastn = str(input("\nUse BLASTN default local scoring system? [y/n]: "))
        if blastn == "y":
            aligner.mode = 'local'
            aligner.match_score = 2
            aligner.mismatch_score = -3
            aligner.open_gap_score = -7
            aligner.extend_gap_score = -2
            print("\nAlignment algorithm:", aligner.algorithm)
            return

        print("\nPlease manually input your scoring preferences.\n")

        mat = str(input("\nUse DNAfull scoring matrix? [y/n]: "))

        if mat == "y":
            aligner.substitution_matrix = dnafull
        else:
            aligner.match_score = float(input("\nCurrent match score: " + 
                str(aligner.match_score) + 
                "\nSet match score: "))
            aligner.mismatch_score = float(input("\nCurrent mismatch score: " + 
                str(aligner.mismatch_score) + 
                "\nSet mismatch score: "))

        aligner.internal_open_gap_score = 0
        aligner.internal_extend_gap_score = 0
        aligner.left_open_gap_score = 0
        aligner.left_extend_gap_score = 0
        aligner.right_open_gap_score = 0
        aligner.right_extend_gap_score = 0

        aligner.internal_open_gap_score = float(
            input("\nCurrent internal open gap score: " + 
            str(aligner.internal_open_gap_score) + 
            "\nSet internal gap open score: "))
        aligner.internal_extend_gap_score = float(
            input("\nCurrent internal extend gap score: " + 
            str(aligner.internal_extend_gap_score) + 
            "\nSet internal gap extend score: "))
        aligner.left_open_gap_score = float(
            input("\nCurrent left open gap score: " + 
            str(aligner.left_open_gap_score) + 
            "\nSet left open gap score: "))
        aligner.left_extend_gap_score = float(
            input("\nCurrent left extend gap score: " + 
            str(aligner.left_extend_gap_score) + 
            "\nSet left extend gap score: "))
        aligner.right_open_gap_score = float(
            input("\nCurrent right open gap score: " + 
            str(aligner.right_open_gap_score) + 
            "\nSet right open gap score: "))
        aligner.right_extend_gap_score = float(
            input("\nCurrent right extend gap score: " + 
            str(aligner.right_extend_gap_score) + 
            "\nSet right extend gap score: "))

        divide = str(input("\nDivide score by length of shorter sequence? [y/n]: "))
        if divide == "y":
            divide_bool = True
        else:
            divide_bool = False
        # divide_bool is global and will be used later with the scoring system

        print("\nAlignment algorithm:", aligner.algorithm)
        
    else:
        aligner.substitution_matrix = dnafull
        aligner.internal_open_gap_score = -10
        aligner.internal_extend_gap_score = -0.5
        aligner.left_open_gap_score = 0
        aligner.left_extend_gap_score = 0
        aligner.right_open_gap_score = 0
        aligner.right_extend_gap_score = 0
        print("\nAlignment algorithm:", aligner.algorithm)
##### Methods

def faToSeqRecordList(input : str) -> list:
    '''
    This function takes in a string representing the fasta filename. It is 
    necessary that the filename is either a complete path or is present in the 
    working directory.
     -> list[SeqRecord]
    '''
    out = []
    for record in SeqIO.parse(input, "fasta"):
        out.append(record)
    return out

def dirToSeqRecordList(input : str) -> list:
    '''
    This function takes in a string representing a directory name. It then
    scans the directory for all sequence files. It reads each and returns
    a list of SeqRecord objects.
    -> list[SeqRecord]
    '''
    output = []
    filenames = glob.glob(input + "/*.seq")
    for filename in filenames:
        with open(filename, 'r') as f:
            data = f.readlines()
        string = ""
        for line in data:
            string = string + (line.upper()).replace("\n", "")
        seq = Seq(string)
        newSeqRec = SeqRecord(seq)
        newSeqRec.id = (str(os.path.basename(filename)).split("_"))[0]
        output.append(newSeqRec)
    return output

def addReverseComplement(input : list) -> list:
    '''
    This function takes in a list of SeqRecord objects. It then appends
    new SeqRecords representing the reverse complement of each SeqRecord to
    the input list.
    '''
    output = []
    for seqrec in input:
        output.append(seqrec)
        newseqrec = SeqRecord((seqrec.seq).reverse_complement())
        newseqrec.id = seqrec.id + "_reverse_complement"
        output.append(newseqrec)
    
    return output

def formatAlignment(alignment : Align.PairwiseAlignment) -> str:
    '''
    This function converts a PairwiseAlignment object into a readable, wrapped
    string. Here, we wrap after 60 characters.
    '''
    lst = format(alignment).split("\n")
    query = lst[0]
    bam = lst[1]
    reference = lst[2]
    length = max(len(bam), len(reference), len(query))
    quotient = math.floor(length / 60)
    output = ""
    for i in range(quotient):
        output += query[0:60] + "\n" + bam[0:60] + "\n" + \
            reference[0:60] + "\n\n"
        query = query[60:]
        bam = bam[60:]
        reference = reference[60:]
    output += query + "\n" + bam + "\n" + reference
    return output

def getAlignment(query : Seq, reference : Seq) -> str:
    '''
    This function returns a well formatted string representing the alignment
    of the two input sequences.
    '''
    alignments = aligner.align(query, reference)
    return(formatAlignment(alignments[0]))

def bestMatch(query_seqrec : SeqRecord, 
    reference_file : list) -> list:
    '''
    Given a sequence and a list of SeqRecords, this function returns the 
    a list with the following information: query sequence id, 
    reference sequence id, query sequence length, reference sequence length,
    and alignment score.
    '''
    global alignment_list
    global divide_bool

    output = []
    output.append(query_seqrec.id)

    best_score = -1
    best_seqrec = SeqRecord(seq=Seq(''))

    query_seqrec_len = len(query_seqrec)
    best_seqrec_len  = 0
    for seqrec in reference_file:
        seqrec_len = len(seqrec)
        length = min(query_seqrec_len, seqrec_len)
        if divide_bool:
            score = aligner.score(query_seqrec.seq, seqrec.seq) / length
        else:
            # print(query_seqrec.seq + "\n" + seqrec.seq)
            score = aligner.score(query_seqrec.seq, seqrec.seq)
            # print(score)

        if score > best_score:
            best_score = score
            best_seqrec = seqrec
            best_seqrec_len = seqrec_len

    
    output.append(best_seqrec.id)
    output.append(query_seqrec_len)
    output.append(best_seqrec_len)
    output.append(best_score)

    name = query_seqrec.id + "_vs_" + best_seqrec.id + ".txt"
    alignment = getAlignment(query_seqrec.seq, best_seqrec.seq)
    alignment_list.append([name, alignment])

    print("\nID: " + output[0] + "\tBest Match: " + output[1] + 
        "\tScore: " + str(output[4]) + "\nTime elapsed:", 
            timedelta(seconds= time.time() - t))

    return output

def calc_df(query_dir : str, reference_file : str):
    '''
    This function calls bestMatch for every SeqRecord present in the
    list generated from the query folder
    '''
    global df
    global check_reverse_complement
    try:
        reference = faToSeqRecordList(reference_file)
    except:
        reference = dirToSeqRecordList(reference_file)

    if check_reverse_complement:
        reference = addReverseComplement(reference)

    try:
        query = faToSeqRecordList(query_dir)
    except:
        query = dirToSeqRecordList(query_dir)
    
    for query_seqrec in query:
        df.loc[len(df)] = bestMatch(query_seqrec, reference)

# Set up the PairwiseAligner object by asking for user input
setUpAligner()

# Get time for calculating elapsed time
t = time.time()

# Get date
d = date.today()

# Compute
calc_df(query_input_path, reference_input_path)

# Make new directory to store results
new_dir_name = path +  "/pC_" + \
    os.path.basename(query_input_path) + \
    "_" + os.path.basename(reference_input_path.split(".")[0]) + "_"+ \
    d.strftime("%m%d%y")

os.mkdir(new_dir_name)
os.chdir(new_dir_name)

# Reset the path to new working directory
path = os.getcwd()

# Add the dataframe as a csv file to the output folder
df.to_csv("summary.csv", sep = "\t", index = False)

# Make new directory to store alignment text files specifically
os.mkdir(path + "/alignments")
os.chdir(path + "/alignments")

# Add text files to new folder
for list in alignment_list:
    with open(list[0], "w") as text_file:
        text_file.write(list[1])
        
# Program completion output notes 
print("\n\n\nDone.\n\n\n")
print("The results were saved to a directory named: \n" + 
    path + "\n")
print("\nTime elapsed:", timedelta(seconds= time.time() - t))
