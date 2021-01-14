import os
import sys
import random
import argparse
from argparse import RawTextHelpFormatter
from itertools import groupby

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class Colors(dict):
	def __init__(self, n=0):
		self.n = n
	def get(self, value):
		if value not in self:
			self[value] = (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1))
		return self[value]

		

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

def get_args():
	usage = 'plot.py [-opt1, [-opt2, ...]] decoded.tsv'
	parser = argparse.ArgumentParser(description='plot', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file in tab separated format')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write the output [stdout]')
	return parser.parse_args()
args = get_args()


fig, ax = plt.subplots()

names = dict()
labels = list()
colors = Colors()
with open(args.infile) as fp:
	for y,line in enumerate(fp):
		for x,name in enumerate(line.rstrip().split('\t')):
			if not x:
				labels.append(name)	
			elif name != '-':
				ax.scatter(x,y,s=100, color=colors.get(name))
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels)	
ax.set_xlabel('gene order')

ax.legend(prop={'size': 6},
		  bbox_to_anchor=(0., 1, 1, 0.1),
		  ncol=len(colors),
		  loc=3,
		  mode="expand",
		  borderaxespad=0.,
		  #handles=[mpatches.Patch(color=col, label=str(lab)) for lab,col in colors.items()])
		  handles=[plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=col, label=str(lab))[0] for lab,col in colors.items()])
fig.savefig('figure.png', bbox_inches='tight')




