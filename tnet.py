#!/usr/bin/python3
# Here is the implementstion of TNet in python3
# Saurav Dhar

# Library imports
from Bio import Phylo
import sys, os

# Global variables




def initialize(input_file):
	print(input_file)
	input_tree = Phylo.read('test/RAxML_rootedTree.9', 'newick')
	print('Is tree bifurcating?', input_tree.is_bifurcating())




def main():
	if not len(sys.argv) == 3:
		raise IndexError("Usage: python3 tnet.py [input phylogeny file] [desired output file]")

	INPUT_TREE_FILE = os.path.abspath(sys.argv[1])
	OUTPUT_FILE = os.path.abspath(sys.argv[2])

	initialize(INPUT_TREE_FILE)


if __name__ == "__main__": main()
