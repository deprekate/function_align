# function_align
Code to multiple align a sets of ordered genes by their function 

The workflow is divided into three steps:
1. encode gene names into letters
2. do multiple alignment with clustalw
3. decode the letters back into gene names


Provided in this repo are several sample phage genomes in GenBank format and a namesfile with several common phage genes 


## Quickstart
```
ls gb/* | xargs -i python3 encode.py {} names.tsv > encoded.fasta
python3 align.py encoded.fasta
python decode.py encoded_aligned.fasta
```
------
## Methods
### 1
For the first step, use the script `encode_genbank.py` to take a genbank file and a namesfile.  The script parses the GenBank file looking for CDS *features* (protein-coding genes).  For each *feature* the script looks for one of the words from the namesfile (with priority given to words at the top of the namesfile) in every *qualifier* field of that *feature*.
To run the script on one genoms you can use the following command:
```
encode.py gb/AB045978.gb names.tsv
```
One of the issues with aligning phage genomes is that many are circular so currently the script looks for an integrase and sets that as the beginning.  The other issue is in which direction (forward or reverse) to start adding to the gene order. Currently the script orders the direction by looking at the order of the collar and terminase genes

### 2
The second step is to do a multiple sequence alighment of all the encoded genenames. This is done using clustalw and the provided custom protein weight matrix, where the weights for matches are 9, mismatches are -9, and wildcards (\*) are 0

The specific clustalw command is 
```
#clustalw2 -INFILE=infile.fasta -OUTFILE=outfile.fasta -OUTPUT=fasta -MATRIX=MATRIX -TYPE=PROTEIN -GAPOPEN=0 -GAPEXT=0 -ALIGN -NEGATIVE
```

### 3
