# Running Environment

The script was written and run only and completely on the symtem: **Windows 7 64 bit professional SP1**. The script is written in Python 3. Since Python is a cross-platform programming language, the script should be able to be run on other platform such Linux and MacOS. If the script or any other files doesn't work on other platforms, perhaps the newline symbols needs to be converted, using NotePad++ or something else.



The program and file directory were "tar and gzip"ed using 7-Zip [https://www.7-zip.org/](https://www.7-zip.org/) on **Windows 7 64 bit professional SP1**. It is easy to "tar and gzip" on Linux server, however, uploading, processing, and downloading are troublesome. it is tricky to "tar and gzip" on Windows system. But the trick is easily resolved here: [https://www.random-host.com/blogs/creating-targz-file-7-zip](https://www.random-host.com/blogs/creating-targz-file-7-zip)



# Dependency

**Python3** needs to be installed. Please refer to [https://www.python.org/](https://www.python.org/)

**BioPython**, a Python package for biologists, needs to be installed. Please refer to [https://biopython.org/](https://biopython.org/)

**BLAST+**, the program to run blast, needs to be installed on the desired platform. The executables of BLAST+ are supposed to automatically added to the PATH environment. Please refer to [http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)



# Execution

The script was written using an IDE called JupyterLab ([https://jupyter.org](https://jupyter.org)). So the raw script file is **main.ipynb**. This file can be opened and run by Jupyter Notebook or JupyterLab. To facilitate non-Jupyter users, a generic Python3 executable script was exported to a file called **main.py**. main.py can be run on different IDEs and platforms.



The input and output paths are hard-coded into the script. If all the input files and the script file are within the same parental directory, the setup will work. This way it is more conveniently to see the results by running each line or cell of codes.



# Choice considerations

In the script, there are many choices to design the code differently. In actual research or production environment, the programmer should understand the whole experimental design and technologies used such as library preparation and sequencing protocol, etc. So the choices can be chosen optimally, dynamically, and interactively.

However, the choices are made based on the best educated guess. For example, should the read sequences first trimmed or first filtered? Should the BLAST+ mask the low complexity subject sequences? Should the primer hits on the reverse strand also be considered? Should the discarded reads also be reported in the additional question 3.



# Files summary

## Input files

test.fna and test.qual

primer and adaptor sequence are hard-coded into the script.

## Result files

Question 1-6: summary.txt

Additional question 1: primer_blast_result.txt and adaptor_blast_result.txt

Additional question 2: filtered_trimmed_test.fna

Additional question 3: read_id_primer_adaptor_hit.txt

## Auxiliary files

All other files.





