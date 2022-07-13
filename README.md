# nema_blast
**nema_blast.py** is a tool for matching sequences based on pairwise alignments. 

Given a list of sequences to be queried and another list of sequences to be referenced, **nema_blast.py** will match each sequence from the query sequence list to a corresponding sequence from the reference sequence list with which it is best aligned.

## Software packages expected in user's PATH
* [Anaconda (includes numpy, pandas)](https://www.anaconda.com/products/distribution)
  * [numpy](https://numpy.org/install/)
  * [pandas](https://pandas.pydata.org/docs/getting_started/install.html)
* [Biopython](https://biopython.org/wiki/Download)
## Installation
    git clone https://github.com/AndersenLab/nema_blast.git
## Usage
**nema_blast.py** analyzes two files/directories and expects to receive two files/directory paths when called.

    python3 nema_blast.py <query_path> <reference_path>

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
    
## Scoring
