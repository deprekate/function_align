import os
import sys
from subprocess import Popen, PIPE, STDOUT

if len(sys.argv) != 2:
	raise ValueError("usage: python3 align.py encoded.fasta")

basename, ext = os.path.splitext(sys.argv[1])
output = b''

try:
	output = Popen(['clustalw','-INFILE=' + sys.argv[1], '-OUTFILE='+basename+'_aligned'+ext, '-OUTPUT=fasta', '-MATRIX=MATRIX', '-TYPE=PROTEIN', '-GAPOPEN=0', '-GAPEXT=0', '-ALIGN', '-NEGATIVE'], stdout=PIPE, stdin=PIPE, stderr=PIPE).stdout.read()
except:
	raise Exception("error: clustalw not found")

with open(basename + '.log', 'w') as fp:
	fp.write(output.decode())
