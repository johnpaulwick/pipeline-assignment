RBIF 100 Assignment 3
Author: John Paul Wick
27 February 2024

This assignment consists of one python script located in /home/wickj/week5, titled pipeline.py.
The purpose of this script is to identify mutations in a series of sequences for molds found on the ears of 50 patients. The molds were seen to be the same species, D. gorgon, but differed in color. This script identifies, for each patient, the mutation present in the sequence reads, which supposedly causes the mold to have its unique color. 

This script reads information from the files present within /home/rbif/week5, and so successful execution relies on its ability to access those files in those locations.
Output is written using relative paths, and so the script can be executed from any directory with output appearing locally. It can be run more than once, without changing or duplicating any of the output. 

This assignment outputs two directories and one text file:
(1) The fastqs directory contains 50 fastq files of the all the sequence data for each patient (each fastq block is one sequence read). Each file is named for its patient.
(2) The bams directory contains a sorted.bam file and a sorted.bam.bai file for each patient (100 files total) from which the pileup is run.
(5) The text file, report.txt, contains a summary of the findings. It includes, for each of the 50 patients: the name of the patient, the color of the mold found on the patient, the number of sequence reads of that mold, the number and percentage of those reads with the mutation, the position of the mutation, the base found at that position in the wildtype sequence, and the base found at that position in the mutated sequences.

Instructions for executing the script are as follows:

(1) Enter an empty directory from which you want to execute the script. Issue the following command to copy in the script:
cp /home/wickj/week5/pipeline.py pipeline.py

(2) Copy and paste the following command to run the script and generate the report:
./pipeline.py

NOTE: The script runtime is less than a minute, but not immediate. Updates at specific progress points in the script's execution are printed in the terminal and a final message will print in the terminal when the final report is ready for viewing.