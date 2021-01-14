import os
import re
import sys
from pathlib import Path
import textwrap
import argparse
from argparse import RawTextHelpFormatter
from itertools import groupby

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def is_valid_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def get_args():
    usage = 'encode_genbank.py [-opt1, [-opt2, ...]] infile'
    parser = argparse.ArgumentParser(description='encode_genbank', formatter_class=RawTextHelpFormatter, usage=usage)
    parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')
    parser.add_argument('namefile', type=is_valid_file, help='names file with one name per line')
    parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
    parser.add_argument('-s', '--squeeze', action="store_true", help='squeeze identical consecutive names into one letter')
    parser.add_argument('-i', '--ignore', action="store_false", help='ignore genes that do not match any of the names')
    parser.add_argument('-e', '--encodings', action="store_true", help='show the names to letter encodings')
    return parser.parse_args()
args = get_args()


names = dict()
with open(args.namefile) as fp:
	for line,letter in zip(fp, 'ARHDCQEGNILKMFPSTWYVBOZX'):
		names[line.rstrip().lower()] = letter

if args.encodings:
	sys.stderr.write(str(names) + "\n")
	exit()

infile = open(args.infile, 'r')
for record in SeqIO.parse(infile, "genbank") :
	letters = []
	for feature in record.features :
		if feature.type=="CDS" :
			letter = '-' if args.ignore else ''
			values = " ".join(map(str,feature.qualifiers.values()))
			for name in names:
				if re.search(rf"\b{name}\b", values, re.IGNORECASE):
					letter = names[name]
					break
			letters.append(letter)

	args.outfile.write(">" + Path(args.infile).stem + '\n')
	# recircularize the order so that integrase is first
	index = letters.index(names['integrase'])
	if letters.index(names['terminase']) < letters.index(names['portal']):
		letters = letters[index:] + letters[:index]
	else:
		letters = list(reversed(letters[:index+1])) + list(reversed(letters[index+1:]))
	if args.squeeze:
		letters = [i[0] for i in groupby(letters)]
	args.outfile.write("".join(letters))
	args.outfile.write('\n')
					

