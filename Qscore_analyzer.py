#!usr/bin/python

'''This program takes a fastq file and graphs the average, median, standard deviation, and variance for all q scores at each nucleotide position'''

import gzip
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description = 'k-merize your fastq files, and generate a graph')
parser.add_argument('-r' , "--read", type = int, nargs = 1, help='read length (base pairs)')
parser.add_argument('-f' , "--file", type = str, nargs = 1, help='Fastq file path')
parser.add_argument('-o' , "--output", type = str, nargs = 1, help='output file (chart) name and/or path')
args = parser.parse_args()

READ_LENGTH = args.read[0]
FILE = args.file[0]
OUTPUT_FILE_NAME = args.output[0]

print('Analyzing the file,', FILE)
print('Read length set to', READ_LENGTH)

READ1_FILE = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz' #This is read1
READ1_INDEX_I7_FILE = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz' #This is i7 index (corresponds to read 1)
READ2_INDEX_I5_FILE = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz' #This is  i5 index (corresponds to read 2)
READ2_FILE = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz' #This is read 2



def init_list(list, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with read length values of 0.0. Must use list = [] first to create empty list. Use list = init_list(list) to run'''
    for i in range(0, READ_LENGTH):
        list.append(value)
    return list

def convert_phred(letter):
    """Converts a single character into a phred score. NOTE: ASSUMES PHRED 33"""
    phred = ord(letter) - 33
    return phred 

def init_list_with_lists(list):
    for x in range(READ_LENGTH):
        zamboni = []
        list.append(zamboni)
    return(list)

def populate_list(file):
    """This function sums together the quality scores (numeric) for each nucleotide, indexed in a fastq file and outputs those
    values in a list of sums"""
    line_counter = 0
    mean_scores = []
    mean_scores = init_list(mean_scores)
    with gzip.open(file) as fh:
        for line in fh:
            line_counter += 1
            if line_counter % 4 == 0:
                letter_counter = 0
                for letter in line.rstrip():
                    # score = convert_phred(letter) #This is necessary when the input file isn't binary (zipped)
                    score = letter -33
                    mean_scores[letter_counter] += score
                    letter_counter += 1
            # if line_counter % 1000000 == 0:
            #     print(str(line_counter) + ' lines complete') #progress bar for short files
    return(mean_scores, line_counter)

mean_scores = []
mean_scores = init_list(mean_scores)

mean_scores, NR = populate_list(FILE)
print(mean_scores)
# print(NR) #prints number of lines in input file

BP_counter = 0
for BP_position in mean_scores:
    mean_scores[BP_counter] = mean_scores[BP_counter] / (NR/4)
    BP_counter += 1


# print(mean_scores)

# LINECOUNT = 0
# print("# Base Pair    Mean Quality Score")
# for x in range(len(mean_scores)):
#     print(LINECOUNT,"\t", mean_scores[x]) 
#     LINECOUNT += 1

BP_POSITION_LIST = []
BP_POSITION_LIST = init_list(BP_POSITION_LIST)
for x in range (READ_LENGTH+1):
    BP_POSITION_LIST[x-1] = x
# print(BP_POSITION_LIST)

plt.plot(BP_POSITION_LIST, mean_scores)
plt.xlabel('Base-pair position')
plt.ylabel('Average Q score')
plt.title(('BP position vs. Average Q score for' + str(FILE)))
plt.savefig(OUTPUT_FILE_NAME)

