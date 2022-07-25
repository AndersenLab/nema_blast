# nema_blast
**nema_blast.py** is a tool for matching nucleotide sequences based on pairwise alignment scores. 

Given a list of sequences to be queried and another list of sequences to be referenced, **nema_blast.py** will match each sequence from the query sequence list to a corresponding sequence from the reference sequence list with which it is best aligned.

**nema_blast.py** may be more efficient at identifying the species associated with an unknown sequence than using [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) if a list of possible sequence matches (reference sequence list) can be provided.

## Software packages expected in user's PATH
* [Anaconda (includes numpy, pandas)](https://www.anaconda.com/products/distribution)
  * [numpy](https://numpy.org/install/)
  * [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
* [Biopython](https://biopython.org/wiki/Download)
## Installation
    git clone https://github.com/AndersenLab/nema_blast.git
## Usage
**nema_blast.py** analyzes two files/directories. As such, **nema_blast.py** expects to receive two files/directory paths when called.

    python3 nema_blast.py <query_path> <reference_path>

An input file must be a properly formatted fasta (.fa) file. An input directory must contain at least one properly formatted sequence (.seq) file.

If the user does not specify a reference file/directory (provides only one file/directory path), then **nema_blast.py** will use the fasta file **reference.fa** as the reference list of sequences.

    python3 nema_blast.py <query_path>
    
Alternatively, the user can choose to call **nema_blast.py** without additional file/directory paths. In this case, **nema_blast.py** will prompt the user to input a query file/directory path and a reference file/directory path.
    
    python3 nema_blast.py

When the program begins, the user will be prompted to provide additional information regarding the alignment scoring system.

First, **nema_blast.py** will ask the user whether it should consider reverse complements of the reference sequences. If the user inputs "y", then **nema_blast.py** will consider, in addition to the existing sequences, the reverse complement of every sequence present in the reference file/directory.
If the strand/direction of all sequences in the reference sequence list are known, the user should not input "y" in order to speed up the program's runtime.

    Consider reverse complements? [y/n]: 
    
Next, **nema_blast.py** will ask the user whether it should use the "EMBOSS Needle Default Scoring System". I recommend using this scoring system for most tasks. The details and source of this scoring system can be found in the next section.

    Use EMBOSS Needle default scoring system? [y/n]: 

If the user inputs "n", then **nema_blast.py** will ask the user whether it should use the "BLASTN default local scoring system". This scoring system implements a local alignment algorithm, which may not be suitable for many tasks. The details and source of this scoring system can be found in the next section.

    Use BLASTN default local scoring system? [y/n]: 
    
If the user inputs "n", then **nema_blast.py** will ask the user to manually input a series of scoring system settings. The next section will provide information regarding the 

    Use DNAfull scoring matrix? [y/n]: 
    # if "n"
    # Set match score: 
    # Set mismatch score: 
    
    Set internal gap open score: 
    Set internal gap extend score: 
    Set left open gap score: 
    Set left extend gap score: 
    Set right open gap score: 
    Set right extend gap score: 
    Divide score by length of shorter sequence? [y/n]: 
    
**nema_blast.py** will create a new directory containing a subdirectory named **alignments**. The subdirectory **alignments** houses formatted individual alignments as well as a file named **summary.csv** which summarizes the best matches and includes additional information regarding sequence length and alignment score.

Here is an example of using **nema_blast.py** to search the **test_data** query directory against the reference file **reference.fa**. 
    
    python3 nema_blast.py test_data bin/reference.fa
    
## Scoring
Please reference the [Biopython documentation](http://biopython.org/DIST/docs/tutorial/Tutorial.html) for information regarding pairwise alignment scoring. It is important to understand how mismatches and gaps are handled by Biopython. Information can be found in section 6.6.2 and section 6.7.

The [DNAfull (also known as EDNAFULL) scoring matrix](https://rosalind.info/glossary/dnafull/) is commonly used for DNA and RNA alignment problems:

$$
\begin{bmatrix}
  &  A  & T &  G &  C &  S &  W  & R &  Y &  K &  M &  B &  V &  H &  D &  N \\
A &  5 & -4 & -4 & -4 & -4 &  1  & 1 & -4 & -4 &  1 & -4 & -1 & -1 & -1 & -2 \\
T & -4 &  5 & -4 & -4 & -4 &  1  &-4 &  1 &  1 & -4 & -1 & -4 & -1 & -1 & -2 \\
G & -4 & -4 &  5 & -4 &  1 & -4  & 1 & -4 &  1 & -4 & -1 & -1 & -4 & -1 & -2 \\
C & -4 & -4 & -4 &  5 &  1 & -4  &-4 &  1 & -4 &  1 & -1 & -1 & -1 & -4 & -2 \\
S & -4 & -4 &  1 &  1 & -1 & -4  &-2 & -2 & -2 & -2 & -1 & -1 & -3 & -3 & -1 \\
W &  1 &  1 & -4 & -4 & -4 & -1  &-2 & -2 & -2 & -2 & -3 & -3 & -1 & -1 & -1 \\
R &  1 & -4 &  1 & -4 & -2 & -2  &-1 & -4 & -2 & -2 & -3 & -1 & -3 & -1 & -1 \\
Y & -4 &  1 & -4 &  1 & -2 & -2  &-4 & -1 & -2 & -2 & -1 & -3 & -1 & -3 & -1 \\
K & -4  & 1 &  1 & -4 & -2 & -2  &-2 & -2 & -1 & -4 & -1 & -3 & -3 & -1 & -1 \\
M &  1 & -4 & -4 &  1 & -2 & -2  &-2 & -2 & -4 & -1 & -3 & -1 & -1 & -3 & -1 \\
B & -4 & -1 & -1 & -1 & -1 & -3  &-3 & -1 & -1 & -3 & -1 & -2 & -2 & -2 & -1 \\
V & -1 & -4 & -1 & -1 & -1 & -3  &-1 & -3 & -3 & -1 & -2 & -1 & -2 & -2 & -1 \\
H & -1 & -1 & -4 & -1 & -3 & -1  &-3 & -1 & -3 & -1 & -2 & -2 & -1 & -2 & -1 \\
D & -1 & -1 & -1 & -4 & -3 & -1  &-1 & -3 & -1 & -3 & -2 & -2 & -2 & -1 & -1 \\
N & -2 & -2 & -2 & -2 & -1 & -1  &-1 & -1 & -1 & -1 & -1 & -1 & -1 & -1 & -1
\end{bmatrix}
$$

Details of the [EMBOSS Needle default scoring system](https://www.ebi.ac.uk/Tools/psa/emboss_needle/):
* scoring matrix: DNAfull
* match: N/A
* mismatch: N/A
* internal open gap: -10
* internal extend gap: -0.5
* right open gap: 0
* right extend gap: 0
* left open gap: 0
* left extend gap: 0
* local/global alignment: global

Details of the [BLASTN default local scoring system](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec118) (as described in the Biopython documentation):
* scoring matrix: N/A
* match: 2
* mismatch: -3
* internal open gap: -7
* internal extend gap: -2
* right open gap: -7
* right extend gap: -2
* left open gap: -7
* left extend gap: -2
* local/global alignment: local
