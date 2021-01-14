import os
import sys
import textwrap
import argparse
from argparse import RawTextHelpFormatter
import tempfile
from subprocess import Popen, PIPE, STDOUT

def is_valid_file(x):
    if not os.path.exists(x):
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def get_args():
    usage = 'align.py [-opt1, [-opt2, ...]] infile'
    parser = argparse.ArgumentParser(description='align', formatter_class=RawTextHelpFormatter, usage=usage)
    parser.add_argument('infile', type=is_valid_file, help='input file in fasta format')
    parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
    parser.add_argument('-t', '--tree', action="store", default=os.devnull, type=argparse.FileType('w'), help='where to write the tree [devnull]')
    return parser.parse_args()
args = get_args()

#basename, ext = os.path.splitext(sys.argv[1])
output = b''

f = tempfile.NamedTemporaryFile(mode='w+b')
t = tempfile.NamedTemporaryFile(mode='w+b')

try:
	output = Popen(['clustalw','-INFILE=' + args.infile, '-OUTFILE='+f.name, '-NEWTREE='+t.name, '-OUTPUT=fasta', '-MATRIX=MATRIX', '-TYPE=PROTEIN', '-GAPOPEN=0', '-GAPEXT=0', '-ALIGN', '-NEGATIVE'], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
except:
	raise Exception(output.decode())

args.outfile.write(f.read().decode())
args.tree.write(t.read().decode())
