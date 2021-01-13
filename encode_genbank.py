import sys
import re

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord




def translate_record(record):
  """Returns a new SeqRecord with translated sequence."""
  return SeqRecord(seq = record.seq.translate(), \
                   id =  record.id + "_translated", \
                   description = record.description)


if len(sys.argv) != 3:
	exit('USAGE: python3 encode_genbank.py infile.gb genenames.tsv')



names = dict()
with open(sys.argv[2]) as fp:
	for line,letter in zip(fp, 'ARHDCQEGNILKMFPSTWYVBOZX'):
		names[line.rstrip().lower()] = letter

#sys.stderr.write(str(names) + "\n")

infile = open(sys.argv[1], 'r')
for record in SeqIO.parse(infile, "genbank") :
	letters = []
	forward = True
	for feature in record.features :
		if feature.type=="CDS" :
			values = " ".join(map(str,feature.qualifiers.values())).lower()
			letter = '*'
			for name in names:
				if re.search(rf"\b{name}\b", values, re.IGNORECASE):
					letter = names[name]
					break
			letters.append(letter)
	print(">", sys.argv[1])
	# recircularize the order so that integrase is first
	index = letters.index(names['integrase'])
	if letters.index('Q') < letters.index('R'):
		letters = letters[index:] + letters[:index]
	else:
		letters = list(reversed(letters[:index+1])) + list(reversed(letters[index+1:]))
	print("".join(letters))
					

#SeqIO.write(proteins, sys.stdout, 'fasta')
