#!/usr/bin/env python
# coding: utf-8

# In[1]:


# detect file existence and delete a file
# run shell command line to call blast within python
import os

# for exit command
import sys

# to calculate arithmatic mean.
import statistics

# biopython is needed for sequence handling.
import Bio

# read and write fasta files
from Bio import SeqIO

# Read the fna and qual files into an SeqRecord iterator provided by the BioPython package.
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

# Parse the blast result in XML format
from Bio.Blast import NCBIXML


# In[2]:


print(Bio.__version__)


# In[3]:


# The answers to first 6 questions will be written to a file called summary.txt in appending mode
# So we need to delete the existing file at the beginning.
if os.path.isfile("./summary.txt"):
    os.remove("./summary.txt")


# In[4]:


# Read the fna and qual files into an SeqRecord iterator provided by the BioPython package.
paired_fasta_qual_iterator = PairedFastaQualIterator(open("test.fna"), open("test.qual"))

# The list of SeqRecord object will be used throughout this script.
paired_fasta_qual_list = list(paired_fasta_qual_iterator)


# In[5]:


# Question 01: Total number of reads in the original dataset
tally_register_01 = len(paired_fasta_qual_list)
print(tally_register_01)

current_output_text = "Total number of reads in the original dataset: " + str(tally_register_01) + "\n"
summary_output_file = open("summary.txt","a+t")
summary_output_file.write(current_output_text)
summary_output_file.close()


# In[6]:


# Question 02: Total number of reads greater than 100 bp in the original dataset
tally_register_02 = 0

for seq_record in paired_fasta_qual_list:
    if len(seq_record.seq) > 100:
        tally_register_02 += 1
print(tally_register_02)

current_output_text = "Total number of reads greater than 100 bp in the original dataset: " + str(tally_register_02) + "\n"
summary_output_file = open("summary.txt","a+t")
summary_output_file.write(current_output_text)
summary_output_file.close()


# In[7]:


# Question 03: Total number of reads with average quality scores greater than 20 in the original dataset
tally_register_03 = 0

for seq_record in paired_fasta_qual_list:
    # letter_annotations stores the information of phred quality scores.
    phred_scores_of_one_sequence = seq_record.letter_annotations["phred_quality"]
    if statistics.mean(phred_scores_of_one_sequence) > 20:
        tally_register_03 += 1
print(tally_register_03)

current_output_text = "Total number of reads with average quality scores greater than 20 in the original dataset: " + str(tally_register_03) + "\n"
summary_output_file = open("summary.txt","a+t")
summary_output_file.write(current_output_text)
summary_output_file.close()


# In[8]:


# Prepare the fasta file of primer sequence
# Write the primer sequence to the hard disk so that blast can take it as input.
primer_sequence = "CGCCGTTTCCCAGTAGGTCTC"
primer_header_line = ">primer"
primer_output_file = open("primer.fna", "wt")
primer_output_file.write(primer_header_line + "\n")
primer_output_file.write(primer_sequence + "\n")
primer_output_file.close()


# In[9]:


# Prepare the fasta file of adaptor sequence
# Write the adaptor sequence to the hard disk so that blast can take it as input.
adaptor_sequence = "ACTGAGTGGGAGGCAAGGCACACAGGGGATAGG"
adaptor_header_line = ">adaptor"
adaptor_output_file = open("adaptor.fna", "wt")
adaptor_output_file.write(adaptor_header_line + "\n")
adaptor_output_file.write(adaptor_sequence + "\n")
adaptor_output_file.close()


# In[10]:


"""
you can run a local copy of BLAST on your own computer.
This allows you to create custom local databases and run unlimited queries
( limited only by your compute power)
Install BLAST+ from NCBI on your computer
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

output format described here
https://www.ncbi.nlm.nih.gov/books/NBK279684/

No hits found? read this webpage: 
You need two more options: -word_size 11 -dust no
https://bioinformatics.stackexchange.com/questions/4226/blastn-no-hits-found

https://www.biostars.org/p/198383/
Since you are using short queries you should try -task blastn-short with your search. 15 nt is getting awfully close to default word size in blast for nucleotides.

Hit on both forward and reverse read? Refer to:
https://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
"""


# In[11]:


# Use python os.system to call command line to run blast locally.

# Make local database. Subject sequences arise from test.fna
blast_makeblastdb_command = "makeblastdb -dbtype nucl -input_type fasta -in test.fna"
os.system(blast_makeblastdb_command)

# Run blast. Query is primer sequence. Subject is the read sequences.
# outfmt 5 is to produce xml output for biopython to parse easily
# outfmt 6 is to produce xml output for human to read easily
blast_blastn_command = "blastn -out primer_blast_result.xml -outfmt 5 -query primer.fna -db test.fna -evalue 0.001 -task blastn-short -dust no"
os.system(blast_blastn_command)
blast_blastn_command = "blastn -out primer_blast_result.txt -outfmt 6 -query primer.fna -db test.fna -evalue 0.001 -task blastn-short -dust no"
os.system(blast_blastn_command)

# Run blast. Query is adaptor sequence. Subject is the read sequences.
blast_blastn_command = "blastn -out adaptor_blast_result.xml -outfmt 5 -query adaptor.fna -db test.fna -evalue 0.001 -task blastn-short -dust no"
os.system(blast_blastn_command)
blast_blastn_command = "blastn -out adaptor_blast_result.txt -outfmt 6 -query adaptor.fna -db test.fna -evalue 0.001 -task blastn-short -dust no"
os.system(blast_blastn_command)


# In[12]:


# Question 04: Total number of reads with primer sequences
# Parse the primer blastn result using modules in biopython
primer_result_handle = open("primer_blast_result.xml")
primer_blast_record = NCBIXML.read(primer_result_handle)
tally_register_04 = len(primer_blast_record.alignments)
primer_result_handle.close()

current_output_text = "Total number of reads with primer sequences: " + str(tally_register_04) + "\n"
summary_output_file = open("summary.txt","a+t")
summary_output_file.write(current_output_text)
summary_output_file.close()


# In[13]:


# Question 05: Total number of reads with adaptor sequences
# Parse the adaptor blastn result using modules in biopython
adaptor_result_handle = open("adaptor_blast_result.xml")
adaptor_blast_record = NCBIXML.read(adaptor_result_handle)
tally_register_05 = len(adaptor_blast_record.alignments)
adaptor_result_handle.close()

current_output_text = "Total number of reads with adaptor sequences: " + str(tally_register_05) + "\n"
summary_output_file = open("summary.txt","a+t")
summary_output_file.write(current_output_text)
summary_output_file.close()


# In[14]:


# Question 06: Total number of reads with both primer and adaptor sequences
primer_hit_reads_IDs = [x.title for x in primer_blast_record.descriptions]
adaptor_hit_reads_IDs = [x.title for x in adaptor_blast_record.descriptions]

# let's find the intersection of the two list:
tally_register_06 = len(list(set(primer_hit_reads_IDs) & set(adaptor_hit_reads_IDs)))

current_output_text = "Total number of reads with both primer and adaptor sequences: " + str(tally_register_06) + "\n"
summary_output_file = open("summary.txt","a+t")
summary_output_file.write(current_output_text)
summary_output_file.close()


# In[15]:


# Additional question 01: generate the following files: Blast output file in the table format.
"""
In legacy BLAST, the tabular output is called -m 8. 
In BLAST+, the tablular output is called -outfmt 6
Please refer to https://github.com/seqan/lambda/wiki/BLAST-Output-Formats
This question has been addressed in the previous blastn part.
"""


# In[16]:


# Addition question 02: generate the following files:
# Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed.
# The question can be split into two steps.
# First step we trim the reads. Second step we filter the reads.


# In[17]:


# To trim the reads, we define a function to get a list that contains a list of hit alignments on a read.
def get_hsps_list(a_seq_record, a_hit_def_list, a_blast_record):
    
    # seq_record.description is the comlete read identifier in the fna file
    # for example: GGBKTAO02CGUTT length=207 xy=0894_3967 region=2 run=R_2010_04_28_09_04_12_
    
    """
    HSP (high-scoring pair) represents region(s) in the hit sequence that contains significant alignment(s) to the query sequence.
    It contains the actual match between your query sequence and a database entry.
    For more information about HSP, please refer to biopython official tutorial documentation.
    Hereafter, multiple HSP, namely HSPs or hsps are used.
    """
    
    seq_id_index = a_hit_def_list.index(a_seq_record.description)
    hsps_list = a_blast_record.alignments[seq_id_index].hsps
    
    if len(hsps_list) > 2:
        sys.exit("There are more than 2 hits on one subject sequence. That doesn't happen in this test data set!")
    
    return hsps_list


# In[18]:


# Input a list of alignment objects, retrieve the slice ranges and return the slice ranges as a sorted list that contains 1 or 2 list element(s)
def get_slice_range(a_hsps_list):
    
    slice_range_list = []
    for a_hsp_obj in a_hsps_list:
        slice_start = a_hsp_obj.sbjct_start
        slice_end = a_hsp_obj.sbjct_end
        
        # print(slice_start, slice_end)
        slice_range_list.append([slice_start, slice_end])
        
    sorted_slice_range_list = sorted(slice_range_list, key = statistics.mean)

    return sorted_slice_range_list


# In[19]:


# Given a list of 1 or 2 slice ranges, trim the read in SeqRecord object, and return the SeqRecord object with its sequence trimmed.
def get_trimmed_seq(a_slice_range_list, a_seq_record):
    
    # First consider only one alignment per subject read sequence
    if len(a_slice_range_list) == 1:
        alignment_01_left = a_slice_range_list[0][0]
        alignment_01_right = a_slice_range_list[0][1]

        """
        We need to discuss two cases of which strand the primer or adaptor hits on:
        (1) +: we keep -->|  |
        (2) -: we keep |  |<--
        """

        # case (1): hit on forward read
        if alignment_01_left < alignment_01_right:
            if (alignment_01_left + alignment_01_right) / 2 <= len(a_seq_record) / 2:
                # print("The only one hit is on the left half of the read sequence")
                trimmed_seq_record = a_seq_record[alignment_01_right: ]
            else:
                # print("The only one hit is on the right half of the read sequence")
                trimmed_seq_record = a_seq_record[: alignment_01_left -1]
                
        # case(2): hit on reverse read
        else:
            if (alignment_01_left + alignment_01_right) / 2 <= len(a_seq_record) / 2:
                # print("The only one hit is on the left half of the read sequence")
                trimmed_seq_record = a_seq_record[alignment_01_left: ]
            else:
                # print("The only one hit is on the right half of the read sequence")
                trimmed_seq_record = a_seq_record[: alignment_01_right -1]
    
        # print(alignment_01_left, alignment_01_right)
    
    
    # Now consider only two alignments per subject read sequence
    if len(a_slice_range_list) == 2:
        alignment_01_left = a_slice_range_list[0][0]
        alignment_01_right = a_slice_range_list[0][1]
        alignment_02_left = a_slice_range_list[1][0]
        alignment_02_right = a_slice_range_list[1][1]

        """
        We need to discuss four cases of which strand the primer or adaptor hits on:
        (1) ++: we keep -->|  |-->
        (2) +-: we keep -->|  |<--
        (3) -+: we keep <--|  |-->
        (4) --: we keep <--|  |<--
        """

        # case (1): forword strand and forward strand
        if alignment_01_left < alignment_01_right and alignment_02_left < alignment_02_right:
            trimmed_seq_record = a_seq_record[alignment_01_right: alignment_02_left -1]

        # case (2): forword strand and reverse strand
        if alignment_01_left < alignment_01_right and alignment_02_left > alignment_02_right:
            trimmed_seq_record = a_seq_record[alignment_01_right: alignment_02_right -1]

        # case (3): reverse strand and forword strand
        if alignment_01_left > alignment_01_right and alignment_02_left < alignment_02_right:
            trimmed_seq_record = a_seq_record[alignment_01_left: alignment_02_left -1]

        # case (4): reverse strand and reverse strand
        if alignment_01_left > alignment_01_right and alignment_02_left > alignment_02_right:
            trimmed_seq_record = a_seq_record[alignment_01_left: alignment_02_right -1]

        # print(alignment_01_left, alignment_01_right, alignment_02_left, alignment_02_right)
    # print("The length of original read sequence is: ", len(a_seq_record))
    # print("The length of trimmed read sequence is: ", len(trimmed_seq_record))
    return trimmed_seq_record


# In[20]:


# First step: get primers and adaptors trimmed.

# first get a list of read IDs (comlete read identifier in the fna file) for all hit with primer or adaptor.
primer_hit_def_list = [x.hit_def for x in primer_blast_record.alignments]
adaptor_hit_def_list = [x.hit_def for x in adaptor_blast_record.alignments]
# print(len(primer_hit_def_list), len(adaptor_hit_def_list))

# set a series of counters to see whether the script runs correctly.
c1 = 0
c2 = 0
c3 = 0
c4 = 0

# Append the read ID, primer hit range, and primer hit range to an accumulation list.
# This list is for the additional question 03 later.
all_hit_range_list = []

# Append the trimmed requence into an accumulation list.
trimmed_paired_fasta_qual_list = []

for seq_record in paired_fasta_qual_list:
    
    # seq_record.description is the comlete read identifier in the fna file
    # for example: GGBKTAO02CGUTT length=207 xy=0894_3967 region=2 run=R_2010_04_28_09_04_12_
    
    condition_primer_exist = seq_record.description in primer_hit_def_list
    condition_adaptor_exist = seq_record.description in adaptor_hit_def_list
    
    """
    To trim a read requence, there are four cases to discuss
    (1) primer exists and adaptor exists
    (2) primer exists and adaptor doesn't exist
    (3) primer doesn't exist and adaptor exists
    (4) primer doesn't exist and adaptor doesn't exist
    """
    
    # case (1): primer exists and adaptor exists
    if condition_primer_exist == True and condition_adaptor_exist == True:

        current_primer_hsps_list = get_hsps_list(seq_record, primer_hit_def_list, primer_blast_record)
        current_primer_slice_range_list = get_slice_range(current_primer_hsps_list)
        
        if len(current_primer_slice_range_list) > 1:
            sys.exit("There are at least 2 primer hits and at least 1 adaptor hit on one subject sequence. That doesn't happen in this test data set!")
        
        current_adaptor_hsps_list = get_hsps_list(seq_record, adaptor_hit_def_list, adaptor_blast_record)
        current_adaptor_slice_range_list = get_slice_range(current_adaptor_hsps_list)
        
        if len(current_adaptor_slice_range_list) > 1:
            sys.exit("There are at least 1 primer hit and at least 2 adaptor hit on one subject sequence. That doesn't happen in this test data set!")
        
        primer_and_adaptor_slice_range_list = current_primer_slice_range_list + current_adaptor_slice_range_list
        sorted_slice_range_list = sorted(primer_and_adaptor_slice_range_list, key = statistics.mean)
        # print(len(sorted_slice_range_list))
        # print(sorted_slice_range_list)
        
        if len(sorted_slice_range_list) > 2:
            sys.exit("There are more than 2 primer or adaptor hits simultaneously on one subject sequence. That doesn't happen in this test data set!")
        
        current_trimmed_seq_record = get_trimmed_seq(sorted_slice_range_list, seq_record)
        
        current_hit_info_list = [seq_record.description, current_primer_slice_range_list, current_adaptor_slice_range_list]
        
        all_hit_range_list.append(current_hit_info_list)
        
        c1 += 1
        
    # case (2): primer exists and adaptor doesn't exist
    if condition_primer_exist == True and condition_adaptor_exist == False:

        current_primer_hsps_list = get_hsps_list(seq_record, primer_hit_def_list, primer_blast_record)
        current_primer_slice_range_list = get_slice_range(current_primer_hsps_list)
        
        sorted_slice_range_list = sorted(current_primer_slice_range_list, key = statistics.mean)
        # print(len(sorted_slice_range_list))
        # print(sorted_slice_range_list)
        
        if len(sorted_slice_range_list) > 2:
            sys.exit("There are more than 2 primer hits on one subject sequence. That doesn't happen in this test data set!")
        
        current_trimmed_seq_record = get_trimmed_seq(sorted_slice_range_list, seq_record)
        
        current_hit_info_list = [seq_record.description, current_primer_slice_range_list, "NA"]
        
        all_hit_range_list.append(current_hit_info_list)
        
        c2 += 1
        
    # case (3): primer doesn't exist and adaptor exists
    if condition_primer_exist == False and condition_adaptor_exist == True:

        current_adaptor_hsps_list = get_hsps_list(seq_record, adaptor_hit_def_list, adaptor_blast_record)
        current_adaptor_slice_range_list = get_slice_range(current_adaptor_hsps_list)
        
        sorted_slice_range_list = sorted(current_adaptor_slice_range_list, key = statistics.mean)
        # print(len(sorted_slice_range_list))
        # print(sorted_slice_range_list)
        
        if len(sorted_slice_range_list) > 1:
            sys.exit("There are more than 1 adaptor hits on one subject sequence. That doesn't happen in this test data set!")
        
        current_trimmed_seq_record = get_trimmed_seq(sorted_slice_range_list, seq_record)
        
        current_hit_info_list = [seq_record.description, "NA", current_adaptor_slice_range_list]
        
        all_hit_range_list.append(current_hit_info_list)
        c3 += 1
    
    # case (4): primer doesn't exist and adaptor doesn't exist
    if condition_primer_exist == False and condition_adaptor_exist == False:
        
        current_trimmed_seq_record = seq_record
        
        current_hit_info_list = [seq_record.description, "NA", "NA"]
        
        all_hit_range_list.append(current_hit_info_list)
        
        c4 += 1
        
    trimmed_paired_fasta_qual_list.append(current_trimmed_seq_record)

print(c1, c2, c3, c4)
print(c1+c2+c3+c4)
print(len(trimmed_paired_fasta_qual_list))
print(len(all_hit_range_list))


# In[21]:


# Second step: filter fasta files that containing reads greater than 100bp, average read quality scores greater than 20
tally_register_filter = 0

# append the filtered requence into an accumulation list.
filtered_trimmed_paired_fasta_qual_list = []

for seq_record in trimmed_paired_fasta_qual_list:
    
    condition_read_length = len(seq_record.seq) > 100
    
    phred_scores_of_one_sequence = seq_record.letter_annotations["phred_quality"]
    condition_quality_scores = statistics.mean(phred_scores_of_one_sequence) > 20
    
    if condition_read_length and condition_quality_scores:
        tally_register_filter +=1
        filtered_trimmed_paired_fasta_qual_list.append(seq_record)

print(tally_register_filter)
print(len(filtered_trimmed_paired_fasta_qual_list))


# In[22]:


# Now we have a list of SeqRecord objects, we'll write them to a FASTA format file.
# Addition question 02: generate the following files is finished now!!!
# Fasta file containing reads greater than 100bp, average read quality scores greater than 20, primers and adaptors trimmed.
SeqIO.write(filtered_trimmed_paired_fasta_qual_list, "filtered_trimmed_test.fna", "fasta")


# In[23]:


# Addition quation 03: Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer and/or adaptor sequences.

"""
To print the addition question 3 to tab-delimited text file prettily, there are four cases to discuss
(1) primer exists and adaptor exists
(2) primer exists and adaptor doesn't exist
(3) primer doesn't exist and adaptor exists
(4) primer doesn't exist and adaptor doesn't exist

To output each sequence record, we need a list, elements of which are 7 strings.
7 strings make up the header lines.
read_id, primer_hit_01_start, primer_hit_01_end, primer_hit_02_start, primer_hit_02_end, adaptor_hit_start, adaptor_hit_end

"""

# set a series of counters to monitor the script running.
c1 = 0
c2 = 0
c3 = 0
c4 = 0

# create an accumulation list to print the information of each sequence record and hits
id_primer_adaptor_output_list = []

for current_hit_range_list in all_hit_range_list:
    
    # case (1): primer exists and adaptor exists
    if type(current_hit_range_list[1]) == list and type(current_hit_range_list[2]) == list:
        
        read_id = current_hit_range_list[0]
        primer_hit_01_start = str(current_hit_range_list[1][0][0])
        primer_hit_01_end = str(current_hit_range_list[1][0][1])
        primer_hit_02_start = "NA"
        primer_hit_02_end = "NA"
        adaptor_hit_start = str(current_hit_range_list[2][0][0])
        adaptor_hit_end = str(current_hit_range_list[2][0][1])
        
        c1 += 1
    
    # case (2): primer exists and adaptor doesn't exist
    if type(current_hit_range_list[1]) == list and type(current_hit_range_list[2]) != list:
        
        # When there is only one primer hit on the subject read sequence.
        if len(current_hit_range_list[1]) == 1:
            
            read_id = current_hit_range_list[0]
            primer_hit_01_start = str(current_hit_range_list[1][0][0])
            primer_hit_01_end = str(current_hit_range_list[1][0][1])
            primer_hit_02_start = "NA"
            primer_hit_02_end = "NA"
            adaptor_hit_start = "NA"
            adaptor_hit_end = "NA"
            
        # When there are only two primer hit on the subject read sequence.
        if len(current_hit_range_list[1]) == 2:
            
            read_id = current_hit_range_list[0]
            primer_hit_01_start = str(current_hit_range_list[1][0][0])
            primer_hit_01_end = str(current_hit_range_list[1][0][1])
            primer_hit_02_start = str(current_hit_range_list[1][1][0])
            primer_hit_02_end = str(current_hit_range_list[1][1][1])
            adaptor_hit_start = "NA"
            adaptor_hit_end = "NA"
            
        c2 += 1
        
    # case (3): primer doesn't exist and adaptor exists
    if type(current_hit_range_list[1]) != list and type(current_hit_range_list[2]) == list:
        
        read_id = current_hit_range_list[0]
        primer_hit_01_start = "NA"
        primer_hit_01_end = "NA"
        primer_hit_02_start = "NA"
        primer_hit_02_end = "NA"
        adaptor_hit_start = str(current_hit_range_list[2][0][0])
        adaptor_hit_end = str(current_hit_range_list[2][0][1])
        
        c3 += 1
        
    # case (4): primer doesn't exist and adaptor doesn't exist
    if type(current_hit_range_list[1]) != list and type(current_hit_range_list[2]) != list:
        
        read_id = current_hit_range_list[0]
        primer_hit_01_start = "NA"
        primer_hit_01_end = "NA"
        primer_hit_02_start = "NA"
        primer_hit_02_end = "NA"
        adaptor_hit_start = "NA"
        adaptor_hit_end = "NA"
        
        c4 += 1
        
    current_string_output = "\t". join([read_id, primer_hit_01_start, primer_hit_01_end, primer_hit_02_start, primer_hit_02_end, adaptor_hit_start, adaptor_hit_end])
    id_primer_adaptor_output_list.append(current_string_output)

print(current_string_output)
print(c1, c2, c3, c4)
print(c1+c2+c3+c4)
print(len(id_primer_adaptor_output_list))


# In[24]:


# Addition quation 03: Tab de-limited text file containing the read identifiers along with the starting and end positions of the primer and/or adaptor sequences. Now the final output is written to the hard disk.

id_primer_adaptor_header_line = "\t". join(["read_id", "primer_hit_01_start", "primer_hit_01_end", "primer_hit_02_start", "primer_hit_02_end", "adaptor_hit_start", "adaptor_hit_end"])

id_primer_adaptor_output_file = open("read_id_primer_adaptor_hit.txt", "wt")

id_primer_adaptor_output_file.write(id_primer_adaptor_header_line + "\n")

for id_primer_adaptor_string in id_primer_adaptor_output_list:
    id_primer_adaptor_output_file.write(id_primer_adaptor_string + "\n")
id_primer_adaptor_output_file.close()

