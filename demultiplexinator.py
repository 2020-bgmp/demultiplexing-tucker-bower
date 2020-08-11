#!/usr/bin/python

#linecount of the fastq inputs: 1452986940
import argparse
import numpy
import gzip
import itertools
import re #regex

parser = argparse.ArgumentParser(description = 'use options -1f for Read 1')
parser.add_argument('-r' , '--readone' , type = str, nargs = 1, help = '.fastq file containing Read 1')
parser.add_argument('-i' , '--indexone' , type = str, nargs = 1, help = '.fastq file containing Read 2 (index i5, forward)')
parser.add_argument('-I' , '--indextwo' , type = str, nargs = 1, help = '.fastq file containing Read 3 (index i7, reverse)')
parser.add_argument('-R' , '--readtwo' , type = str, nargs = 1, help = '.fastq file containing Read 4')
parser.add_argument('-l' , '--indexlist' , type = str, nargs = 1, help = 'File containing forward (i5) index names and sequences, tab separated')
parser.add_argument('-o' , '--output' , type = str, nargs = 1, help = 'Output directory')
args = parser.parse_args()

R1_FASTQ = gzip.open(args.readone[0], "r") #T
R2_FASTQ = gzip.open(args.indexone[0], "r")
R3_FASTQ = gzip.open(args.indextwo[0], "r")
R4_FASTQ = gzip.open(args.readtwo[0], "r")
INDEX_FILE = open(args.indexlist[0], "r")
OUTPUT_DIR = args.output[0]

def convert_phred(letter):
    """Converts a single character into a phred score. NOTE: ASSUMES PHRED 33"""
    phred = ord(letter) - 33
    return phred 

def index_dictionary_inator():
    '''Creates a dictionary where the keys are the names of the indexes and the values are the index sequences. Dependent on input from
    the -l specified file, which has been opened as INDEX_FILE. Index file must be formatted correctly. Example:
    B1    GATTACAC'''
    index_dict = {}
    while True:
        index = INDEX_FILE.readline()
        if index == '':
            break
        index_name = re.findall('(\w+)\s', index)
        index_sequence = re.findall('\w+\s+(\w+)', index)
        index_dict[index_name[0]] = index_sequence[0]
    return index_dict

INDEX_DICT = index_dictionary_inator()

def reverse_complement_inator(dna_seq):
    '''Given a DNA sequence, returns the reverse complement of that sequence. Acceptable characters, A,T,G,C,N'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(dna_seq)) 
    #.join method joins all the items in a tuple into a string, whatever is in front of .join is the separator
    #.get returns values for a key, the second arg specifies what to return if the key is not in the dict. In this case i just won't translate it
    #In summary, the for loop generates a tuple where each element in the tuple is the complement of the base in the reversed seq. And .join combines the elements of that tuple into a string
    return reverse_complement

REVERSE_INDEX_DICT = INDEX_DICT.copy()
for key in REVERSE_INDEX_DICT:
    REVERSE_INDEX_DICT[key] = reverse_complement_inator(REVERSE_INDEX_DICT[key])
#This chunk of code creates a dict where the keys are the index names and the values are the reverse complement sequences (i7 indexes)
#dependent on INDEX_DICT and reverse_complement_inator execution upstream.

INDEX_PERMUTATIONS_LIST = list(itertools.permutations(INDEX_DICT.keys(), 2))
#In this list, each possible permutation of index names is a tuple containing the two index names.

INDEX_PERMU_COUNTER_DICT = {}
for permu in INDEX_PERMUTATIONS_LIST:
    INDEX_PERMU_COUNTER_DICT[permu] = 0
#In this dict, keys = permutations of index names, values = 0 (will be counted as we iterate the file)
    
LOW_QUALITY_R1_FILE = open(''.join([OUTPUT_DIR, 'low_quality_R1.fastq']), 'w+')
LOW_QUALITY_R4_FILE = open(''.join([OUTPUT_DIR, 'low_quality_R4.fastq']), 'w+')
INDEX_HOPPED_R1_FILE = open(''.join([OUTPUT_DIR, 'index_hopped_R1.fastq']), 'w+')
INDEX_HOPPED_R4_FILE = open(''.join([OUTPUT_DIR, 'index_hopped_R4.fastq']), 'w+')
STATS_FILE = open(''.join([OUTPUT_DIR, 'stats.txt']), 'w+')
INDEX_PERMUTATIONS_FILE = open(''.join([OUTPUT_DIR, 'index_permutations_counts.txt']), 'w+')
#creating output files in output directory

#For the multiplexinator, we need the index dict keys to be the sequences, because keys are hashable and values are not. So we're going to create two new dictionaries,
#from index_dict and reverse_index_dict, but with the keys and values swapped
SWAPPED_INDEX_DICT = dict([(value, key) for key, value in INDEX_DICT.items()]) #swapping keys and values. This would not work if there were multiple values/key
SWAPPED_REVERSE_INDEX_DICT = dict([(value, key) for key, value in REVERSE_INDEX_DICT.items()])

def core_multiplexinator():
    '''This is it. The engine of this script. Sorts all 4 input fastq files into the appropriate 52 output files and counts occurrences of matched indexes, hopped indexes, and low
    quality reads that get thrown out'''
    matched_counter = 0
    hopped_counter = 0
    low_quality_counter = 0
    reads_processed_counter = 0
    #setting up some counters for the stats.txt file
    while True:
        reads_processed_counter += 1
        if reads_processed_counter % 15000000 == 0:
            print(str(reads_processed_counter * 100 / 1500000000) + '%') #Note, this is the only 'hardcode' in here. Need to do a linecount of the file to recreate this progress bar
        R1_header = R1_FASTQ.readline().rstrip().decode() #reads the line, strips the new line character, and gets rid of the weird b in front of the string.
        if R1_header == '':
            break
            #EOF
        R1_sequence = R1_FASTQ.readline().rstrip().decode()
        R1_sep = R1_FASTQ.readline().rstrip().decode()
        R1_qscores = R1_FASTQ.readline().rstrip().decode()
        R2_header = R2_FASTQ.readline().rstrip().decode()
        R2_sequence = R2_FASTQ.readline().rstrip().decode()
        R2_sep = R2_FASTQ.readline().rstrip().decode()
        R2_qscores = R2_FASTQ.readline().rstrip().decode()
        R3_header = R3_FASTQ.readline().rstrip().decode()
        R3_sequence = R3_FASTQ.readline().rstrip().decode()
        R3_sep = R3_FASTQ.readline().rstrip().decode()
        R3_qscores = R3_FASTQ.readline().rstrip().decode()
        R4_header = R4_FASTQ.readline().rstrip().decode()
        R4_sequence = R4_FASTQ.readline().rstrip().decode()
        R4_sep = R4_FASTQ.readline().rstrip().decode()
        R4_qscores = R4_FASTQ.readline().rstrip().decode()
        #grabbing each line of each fastq input file for one read
        R1_output_header = ''.join(R1_header + ' ' + R2_sequence + '-' + R3_sequence)
        R4_output_header = ''.join(R4_header + ' ' + R2_sequence + '-' + R3_sequence)
        #adding the index pair sequences to the output headers
        low_qscore_counter_R2 = 0 
        low_qscore_counter_R3 = 0 
        for score in R2_qscores:
            if convert_phred(score) <= 20:
                low_qscore_counter_R2 += 1
        for score in R3_qscores:
            if convert_phred(score) <= 20:
                low_qscore_counter_R3 += 1
        #Creating a counter of how many qscores below 20 are in the qscores of the indexes. The cutoff in this program is the index reads must not have more than 1 q score below 20
        if low_qscore_counter_R2 >= 3 or low_qscore_counter_R3 >= 3 or R2_sequence not in list(SWAPPED_INDEX_DICT.keys()) or R3_sequence not in list(SWAPPED_REVERSE_INDEX_DICT.keys()):
            #HUGE NOTE HERE. In the index file of indexes I was given, the indexes line up with R3, not R2. I assumed R2 was forward like a fool for so long. Always check to 
            #make sure you know which direction your data is in
            low_quality_counter += 1
            LOW_QUALITY_R1_FILE.write(R1_output_header + '\n')
            LOW_QUALITY_R1_FILE.write(R1_sequence + '\n')
            LOW_QUALITY_R1_FILE.write(R1_sep + '\n')
            LOW_QUALITY_R1_FILE.write(R1_qscores + '\n')
            LOW_QUALITY_R4_FILE.write(R4_output_header + '\n')
            LOW_QUALITY_R4_FILE.write(R4_sequence + '\n')
            LOW_QUALITY_R4_FILE.write(R4_sep + '\n')
            LOW_QUALITY_R4_FILE.write(R4_qscores + '\n')
        #Checks for any of multiple conditions that qualify a read as a low_quality read, if they apply, writes those reads to their appropriate files
        elif SWAPPED_INDEX_DICT[R2_sequence] != SWAPPED_REVERSE_INDEX_DICT[R3_sequence]:
            hopped_counter += 1
            #I think i'm gonna want to make this dict as I go if this doesn't work
            INDEX_PERMU_COUNTER_DICT[(SWAPPED_INDEX_DICT[R2_sequence], SWAPPED_REVERSE_INDEX_DICT[R3_sequence])] += 1
            INDEX_HOPPED_R1_FILE.write(R1_output_header + '\n')
            INDEX_HOPPED_R1_FILE.write(R1_sequence + '\n')
            INDEX_HOPPED_R1_FILE.write(R1_sep + '\n')
            INDEX_HOPPED_R1_FILE.write(R1_qscores + '\n')
            INDEX_HOPPED_R4_FILE.write(R4_output_header + '\n')
            INDEX_HOPPED_R4_FILE.write(R4_sequence + '\n')
            INDEX_HOPPED_R4_FILE.write(R4_sep + '\n')
            INDEX_HOPPED_R4_FILE.write(R4_qscores + '\n')
        #Checks if indexes don't match. if they don't they are sorted into the index hopped fastqs
        elif SWAPPED_INDEX_DICT[R2_sequence] == SWAPPED_REVERSE_INDEX_DICT[R3_sequence]:
            matched_counter += 1
            r1_output = open(''.join([OUTPUT_DIR, 'matched_', SWAPPED_INDEX_DICT[R2_sequence], '_R1.fastq']), 'a')
            r4_output = open(''.join([OUTPUT_DIR, 'matched_', SWAPPED_REVERSE_INDEX_DICT[R3_sequence], '_R4.fastq']), 'a')
            #temporary variables set to open the matched output file. these will be different in every loop depending on the indexes
            r1_output.write(R1_output_header + '\n')
            r1_output.write(R1_sequence + '\n')
            r1_output.write(R1_sep + '\n')
            r1_output.write(R1_qscores + '\n')
            r4_output.write(R4_output_header + '\n')
            r4_output.write(R2_sequence + '\n')
            r4_output.write(R2_sep + '\n')
            r4_output.write(R2_qscores + '\n')
    STATS_FILE.write('Properly matched indexes: ' + str(matched_counter) + '\n')
    STATS_FILE.write('Low quality indexes: ' + str(low_quality_counter) + '\n')
    STATS_FILE.write('Hopped indexes: ' + str(hopped_counter) + '\n')
    STATS_FILE.write('Total records processed: ' + str(reads_processed_counter) + '\n')
    INDEX_PERMUTATIONS_FILE.write(str(INDEX_PERMU_COUNTER_DICT))



core_multiplexinator()
