#!/usr/bin/python

READ1_FILE = 'projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz' #This is read1
READ1_INDEX_I7_FILE = 'projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz' #This is i7 index (corresponds to read 1)
READ2_INDEX_I5_FILE = 'projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz' #This is  i5 index (corresponds to read 2)
READ2_FILE = 'projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz' #This is read 2

import numpy
import gzip

B1 = 'GTAGCGTA' 
A5 = 'CGATCGAT'
C1 = 'GATCAAGG'
B9 = 'AACAGCGA'    
C9 = 'TAGCCATG'    
C3 = 'CGGTAATC'
B3 = 'CTCTGGAT'    
C4 = 'TACCGGAT'    
A11 = 'CTAGCTCA'
C7 = 'CACTTCAC'    
B2 = 'GCTACTCT'    
A1 = 'ACGATCAG'
B7 = 'TATGGCAC'    
A3 = 'TGTTCCGT'    
B4 = 'GTCCTAAG'
A12 = 'TCGACAAG'    
C10 = 'TCTTCGAC'    
A2 = 'ATCATGCG'
C2 = 'ATCGTGGT'    
A10 = 'TCGAGAGT'    
B8 = 'TCGGATTC'
A7 = 'GATCTTGC'    
B10 = 'AGAGTCCA'    
A8 = 'AGGATAGC'

Suquence rever
#with gzip.open(READ1_FILE) as fh:
    