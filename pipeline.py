#!/usr/bin/env python3


# RBIF 100 Assignment 3
# Author: John Paul Wick
# 27 February 2024


# Modules required for pipeline.py:
import gzip
import os
import subprocess
import pysam
import csv


##############################################################################################################################################################


# STEP 1 - PATIENT FASTQ FILES


# FASTQ parser class (copied in from rbif/week5 directory):
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
        
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)


# Creating fastqs directory:
try:
    os.mkdir('fastqs')
except:
    True


# Reading in clinical data:
with open('/home/rbif/week5/harrington_clinical_data.txt') as f:
    reader = csv.reader(f, delimiter="\t")
    patientlist = list(reader)
patientlist.pop(0) # Removes column headers


# Defining arrays:
names = [item[0] for item in patientlist]
colors = [item[1] for item in patientlist]
barcodes = [item[2] for item in patientlist]


# Creating patient files:
for name in names:
    with open(f'fastqs/{name}_trimmed.fastq', 'w') as fp:
        pass


# The function 'writetofileandtrim' adds each fastq block from the pooled data to its corresponding patient's file. It takes two arguments: name and barcode of patient.
# This function also trims barcodes and ends from the sequence and quality score lines. 
# A confirmation message is printed in the terminal upon completion of each patient's fastq file.
def writetofileandtrim(name, barcode):
    fastqfile = ParseFastQ("/home/rbif/week5/hawkins_pooled_sequences.fastq")
    file = open(f'fastqs/{name}_trimmed.fastq', 'w')
    for fastq_obj in fastqfile:
        if fastq_obj[1].startswith(barcode): # This will filter it so only the correct reads are processed and added to the patient's file.
            # Naming each line from each fastq block that belongs to the patient:
            header = fastq_obj[0]
            sequence = fastq_obj[1]
            line3 = fastq_obj[2]
            qualityscores = fastq_obj[3]
            ########################################
            # This next block of code determines where to trim the end off:
            case1 = "DD"
            case2 = "DF"
            case3 = "FD"
            case4 = "FF" 
            pos1 = qualityscores.find(case1)
            pos2 = qualityscores.find(case2)
            pos3 = qualityscores.find(case3)
            pos4 = qualityscores.find(case4)
            positions = [pos1, pos2, pos3, pos4]
            positions = [pos for pos in positions if pos > 0] # This removes negative values from the list (if there is no match in the line, the find function returns -1)
            trimsite = min(positions)
            #########################################
            # Finally, the fastq block is trimmed and added to its patient file:
            file.writelines(f'\
{header}\n\
{sequence[5:trimsite]}\n\
{line3}\n\
{qualityscores[5:trimsite]}\n')
    
    file.close()
    print(f"{name}'s fastq file is written and trimmed!")


# Now I execute the function on each patient:
for patient in patientlist:
    writetofileandtrim(patient[0],patient[2]) # The name and barcode arguments are passed as arguments into the function for each patient


##############################################################################################################################################################


# STEP 2 - BAM FILES


# Creating the bams directory:
try:
    os.mkdir('bams')
except:
    True


# Creating an index of the reference sequence and placing index files inside a new index directory:
try:
    os.mkdir('index')
except:
    True
subprocess.run(['cp', '/home/rbif/week5/dgorgon_reference.fa', 'index/ref.fa'])
subprocess.run(['bwa', 'index', 'index/ref.fa'])


# Aligning each fastq file to the reference and creating a corresponding sam file in the bams directory:
for name in names:
    fastqfile = str(f'fastqs/{name}_trimmed.fastq')
    samfile = str(f'bams/{name}.sam')
    createsams = "bwa mem index/ref.fa " + fastqfile + "> " + samfile
    os.system(createsams)


##############################################################################################################################################################


# STEP 3 - INDEXED SORTED BAM FILES


# Then, I convert each sam file to a bam file, in the same directory:
for name in names:
    samfile = str(f'bams/{name}.sam')
    bamfile = str(f'bams/{name}.bam')
    converttobams = "samtools view -bS " + samfile + "> " + bamfile
    os.system(converttobams)


# Now I will sort the bam files into sorted.bam files:
for name in names:
    bamfile = str(f'bams/{name}.bam')
    sortedbamfile = str(f'bams/{name}.sorted.bam')
    sortbams = "samtools sort -m 100M -o " + sortedbamfile + " " + bamfile
    os.system(sortbams)
    

# Then, I remove the old unneeded sam files and unsorted bam files:
for name in names:
    samfile = str(f'bams/{name}.sam')
    bamfile = str(f'bams/{name}.bam')
    deletesamsandbams = "rm " + samfile + " " + bamfile
    os.system(deletesamsandbams)


# Indexing the sorted bam files into sorted.bam.bai files (with confirmation messages):
for name in names:
    sortedbamfile = str(f'bams/{name}.sorted.bam')
    indexsortedbams = "samtools index " + sortedbamfile 
    os.system(indexsortedbams)
    print(f"{name}'s sorted bam file is complete and ready for pileup!")


##############################################################################################################################################################


# STEP 4 - PILEUP


# I will create an intermediate file with the pileup data for each patient, before generating the final report document.
reportfile = open("pileupdata.txt", "w")
# Adding a header so I know what each column contains:
reportfile.writelines("Name\tColor\tNumber of Reads\tNumber of Reads With Mutation\tPercent of Reads with Mutation\tMutation Position\tWild-type Base\tMutated Base\n")
reportfile.close()


# The following function performs the pileup of a given patient's sorted bam file, and outputs all the necessary information about the mutations from the bam file to a new txt file.
def pileup(name, color):
    samfile = pysam.AlignmentFile(f"bams/{name}.sorted.bam", "rb")
    for pileupcolumn in samfile.pileup():
        # I use a dictionary to count up the bases at each position:
        ntdict = {}
        Acount = 0
        Ccount = 0
        Tcount = 0
        Gcount = 0
        numberofreads = pileupcolumn.n
        position = pileupcolumn.pos + 1 # The nucleotide position is equal to the pythonic index position plus one, because the index positions start at zero.
        for pileupread in pileupcolumn.pileups: # This inner loop is for each individual base in the column currently being iterated.
            if not pileupread.is_del and not pileupread.is_refskip:
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                # I count all the bases in the column:
                if base == "A":
                    Acount = Acount + 1
                elif base == "C":
                    Ccount = Ccount + 1
                elif base == "T":
                    Tcount = Tcount + 1
                elif base == "G":
                    Gcount = Gcount + 1
        # I add the counts for each base in the column to ntdict:
        ntdict["A"] = Acount
        ntdict["C"] = Ccount
        ntdict["T"] = Tcount
        ntdict["G"] = Gcount
        # Now I calculate the percentages of each base for the position's pileup column:
        Apercentage = ntdict["A"] / numberofreads * 100
        Cpercentage = ntdict["C"] / numberofreads * 100
        Tpercentage = ntdict["T"] / numberofreads * 100
        Gpercentage = ntdict["G"] / numberofreads * 100
        # If any of the percentages are NOT 0 or 100, I add the information to a dictionary (the key is the base and the value is the percentage):
        flaggedbasesandpercents = {}
        if Apercentage != 0 and Apercentage != 100:
            flaggedbasesandpercents["A"] = Apercentage
        if Cpercentage != 0 and Cpercentage != 100:
            flaggedbasesandpercents["C"] = Cpercentage
        if Tpercentage != 0 and Tpercentage != 100:
            flaggedbasesandpercents["T"] = Tpercentage
        if Gpercentage != 0 and Gpercentage != 100:
            flaggedbasesandpercents["G"] = Gpercentage
        # If this dictionary is not empty, it means that there are mutations in the column. These columns will be processed, while columns whose dictionaries are empty will be skipped.
        if flaggedbasesandpercents: # The following blocks of code are only run for iterations of columns containing mutations. OTherwise the loop skips to the next iteration.
            referencefasta = open("index/ref.fa", "r")
            referencelines = referencefasta.readlines()
            referencesequence = referencelines[1]
            for key in flaggedbasesandpercents:
                if key != referencesequence[pileupcolumn.pos]:  # It compares the problematic bases in the dictionary to the base in the reference at the position corresponding to the column.
                    mutation = key  # The base which is different is determined to be the mutation.
                    percent = round(flaggedbasesandpercents[key], 2) # The rounded percentage of the base is also named.
                    numberofmutatedreads = ntdict[f"{mutation}"]
                    reference = referencesequence[pileupcolumn.pos] # The wild-type base in the reference strain is also named.
            referencefasta.close()
            # Finally, the output of this function is a string for each column with a mutation containing all relevant information about the mutation and patient.
            # This string is added to the pileupdatafile I created, from which I'll draw information to generate the final report.
            patientdata = f"{name}\t{color}\t{numberofreads}\t{numberofmutatedreads}\t{percent}\t{position}\t{reference}\t{mutation}\n"
            file = open("pileupdata.txt", "a")
            file.writelines(patientdata)
            file.close()

    samfile.close()


# Now I run the pileup on a loop which will, for each patient, add a line corresponding to each position containing mutations in at least some of their sequence reads.
# The file will contain all the information necessary to generate a final report.
for patient in patientlist:
    pileup(patient[0], patient[1])
    print(f"{patient[0]}'s pileup data is complete!") # Progress is documented in the terminal again with confirmation messages.


##############################################################################################################################################################


# STEP 5 - FINAL REPORT


# Creating the text file:
reportfile = open("report.txt", "w")
# Adding a header:
reportfile.writelines("FINAL REPORT: MOLD MUTATIONS BY PATIENT\n\n")
reportfile.close()


# Reading in the data from the temporary pileupdata text file:
with open('pileupdata.txt') as f:
    reader = csv.reader(f, delimiter="\t")
    molddatalist = list(reader)
molddatalist.pop(0) # This removes the headers from the data.


# Defining arrays (these are not used by the script but are present for enhanced readability):
patientnames = [item[0] for item in molddatalist]
moldcolors = [item[1] for item in molddatalist]
reads = [item[2] for item in molddatalist]
mutatedreads = [item[3] for item in molddatalist]
percents = [item[4] for item in molddatalist]
mutationpositions = [item[5] for item in molddatalist]
wildtypebases = [item[6] for item in molddatalist]
mutatedbases = [item[7] for item in molddatalist]


# The following function delivers, for a given patient (as present in molddatalist from the pileupdata), the mold mutation data in a readable format and adds it to the report file.
def writetofinalreport(patient):
    reportfile = open("report.txt", "a")
    reportfile.writelines(f"\
Patient name: {patient[0]}.\n\
The mold sample found on {patient[0]} was the color {(patient[1]).lower()}.\n\
There were {patient[2]} reads of the mold sequence used, with {patient[3]} or {patient[4]}% of the reads containing the mutation.\n\
The mutation in these samples is found at position {patient[5]}, in which the wild-type base {patient[6]} is replaced with {patient[7]}.\n\n")
    reportfile.close()


# On a loop for each patient, I add all the mold mutation data to the report:
for patient in molddatalist:
    writetofinalreport(patient)
reportfile = open("report.txt", "a")
reportfile.writelines("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
reportfile.close()


# I no longer need the index directory and the pileupdata file, so I'll get rid of them:
removeindexdirectory = "rm -r index"
os.system(removeindexdirectory)
removepileupdatafile = "rm pileupdata.txt"
os.system(removepileupdatafile)


# Finally, I print a message in the terminal to alert the user that the report is ready:
print("\
\n\n\n\n\n\n\t\t\t#####################################\n\n\n\t\t\
Process is complete! See report.txt for final report.\
\n\n\n\t\t\t#####################################\n\n\n")


##############################################################################################################################################################