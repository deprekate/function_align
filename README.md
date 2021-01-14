# function_align
Code to multiple-align sets of ordered genes by their function (or words) 

The workflow is divided into three steps:
1. encode gene names into letters
2. do multiple alignment with clustalw
3. decode the letters back into gene names


Provided in this repo are several sample phage genomes in GenBank format and a namesfile with several common phage genes 


## Quickstart
```
ls gb/* | xargs -i python3 encode.py -i {} names.tsv > encoded.fasta
python3 align.py encoded.fasta > aligned.fasta
python decode.py aligned.fasta names.tsv > decoded.tsv
```
------
## Methods
### 1 encode
For the first step, use the script `encode_genbank.py` to take a genbank file and a namesfile.  The script parses the GenBank file looking for CDS *features* (protein-coding genes).  For each *feature* the script looks for one of the words from the namesfile (with priority given to words at the top of the namesfile) in every *qualifier* field of that *feature*.
To run the script on one genoms you can use the following command:
```
python3 encode.py -i gb/AB045978.gb names.tsv
```
The `-i` flag tells the script to ignore genes that are not in the name list

One of the issues with aligning phage genomes is that many are circular so currently the script looks for an integrase and sets that as the beginning.  The other issue is in which direction (forward or reverse) to start adding to the gene order. Currently the script orders the direction by looking at the order of the collar and terminase genes

### 2 align
The second step is to do a multiple sequence alighment of all the encoded genenames. This is done using clustalw and the provided custom protein weight matrix, where the weights for matches are 9, mismatches are -9, wildcards (\*) are 0, and gaps are 0
```
python3 align.py encoded.fasta
```

### 3 decode
The last step is to take the multiuple sequence alignments and decode them back into gene names
```
python3 decode.py aligned.fasta names.tsv
```

### 4 plot
We have supplied a script to plot the tab separated gene name file using matplotlib
```
python3 plot.py decoded.tsv
```

and the figure for the provided example genomes should look like:
![](https://github.com/deprekate/function_align/blob/main/figure.png)
