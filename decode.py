import os
import re
import sys
from pathlib import Path
import textwrap
import argparse
from argparse import RawTextHelpFormatter
from itertools import groupby


def is_valid_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def get_args():
    usage = 'encode_genbank.py [-opt1, [-opt2, ...]] infile namefile'
    parser = argparse.ArgumentParser(description='encode_genbank', formatter_class=RawTextHelpFormatter, usage=usage)
    parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')
    parser.add_argument('namefile', type=is_valid_file, help='names file with one name per line')
    parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
    return parser.parse_args()
args = get_args()


names = dict()
with open(args.namefile) as fp:
	for line,letter in zip(fp, 'ARHDCQEGNILKMFPSTWYVBOZX'):
		names[letter] = line.rstrip().lower()
names['-'] = '-'
#sys.stderr.write(str(names) + "\n")

with open(args.infile, 'r') as fp:
	for line in fp:
		if line.startswith(">"):
			args.outfile.write(line.replace(">", "").rstrip())
		else:
			for letter in line.rstrip():
				args.outfile.write('\t')
				args.outfile.write(names[letter])
			args.outfile.write('\n')

