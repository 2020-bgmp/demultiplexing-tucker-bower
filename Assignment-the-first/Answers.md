cd# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |
| 1294_S1_L008_R2_001.fastq.gz |  |
| 1294_S1_L008_R3_001.fastq.gz |  |
| 1294_S1_L008_R4_001.fastq.gz |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. ![R1 Qscore distribution](https://github.com/2020-bgmp/demultiplexing-tucker-bower/blob/master/Assignment-the-first/Qscore_analyzer_outputs/Mean_Qscore_by_position_R1.png)
    3. ![R2 Qscore distribution](https://github.com/2020-bgmp/demultiplexing-tucker-bower/blob/master/Assignment-the-first/Qscore_analyzer_outputs/Mean_Qscore_by_position_R2.png)
    ![R3 Qscore distribution](https://github.com/2020-bgmp/demultiplexing-tucker-bower/blob/master/Assignment-the-first/Qscore_analyzer_outputs/Mean_Qscore_by_position_R3.png)
    ![R4 Qscore distribution](https://github.com/2020-bgmp/demultiplexing-tucker-bower/blob/master/Assignment-the-first/Qscore_analyzer_outputs/Mean_Qscore_by_position_R4.png)
    
## Part 2
1. Define the problem
Here we have multiplexed data that comes from a parallel run of DNA from many different samples. We need to use the corresponding indexes to sort the samples into their own FASTQ files. Additionally, some indexes are too low quality to be considered in downstream analysis, and will be sorted separately. 

Furthermore, index hopping is when a library DNA molecule has 2 mismatched indexes on either side of the molecule. We would like to analyze index hopping in our data. To address these we will check the indexes when sorting the data into separate fastq files and keep a running count of how many occurrences exist of every possible permutation of indexes. Index hopped reads will also be sorted separately from matched paired ends.
2. Describe output
Outputs:
    -24 output files of matched READ1 indexed reads (one for each index)
    -24 output files of matched READ2 indexed reads
    -2 files of reads that were index hopped (one for each READ)
    -2 files for non-matching or low quality index pairs
    -1 file containing the following: 
        -A count of the number of read-pairs with properly matched indexes
        -A count of the read-pairs with index hopping issues
        -A count of the read-pairs with unknown indexes (q_score is too low, or an N is present)
        -Counts for how many times each possible permutation of indexes happened (ex: B1_B9: 2700)
         Indicates index hopped pair of B1 index on read 1 and B9 on read 2 happened 2700 times
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
