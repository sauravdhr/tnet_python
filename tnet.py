#!/usr/bin/python3
# Here is the implementstion of TNet in python3
# Saurav Dhar

# Library imports
from Bio import Phylo
import sys, os

# Global variables
score = {}
solution_count = {}
hosts = []
transmission_edges = []




def initialize_tree(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True
	print('Terminals:', len(input_tree.get_terminals()))
	print('Nonterminals:', len(input_tree.get_nonterminals()))
	if not input_tree.is_bifurcating():
		raise IndexError("Input tree is not bifurcating.")
	# Phylo.draw_ascii(rooted_tree)
	return input_tree

def initialize_leaf_nodes(rooted_tree):
	temp_host = []
	for terminal in rooted_tree.get_terminals():
		terminal.name = terminal.name.split('_')[0]
		temp_host.append(terminal.name)

	hosts = list(set(temp_host))
	print('Total hosts: ', len(hosts))

	for terminal in rooted_tree.get_terminals():
		temp = []
		count = []
		for host in hosts:
			if host==terminal.name:
				temp.append(0)
				count.append(1)
			else:
				temp.append(9999999999)
				count.append(0)

		score[terminal] = temp
		solution_count[terminal] = count

def initialize_score_count(node):
	l_score = score[node.clades[0]].copy()
	r_score = score[node.clades[1]].copy()
	l_count = solution_count[node.clades[0]].copy()
	r_count = solution_count[node.clades[1]].copy()

	length = len(l_score)
	temp_score = []
	temp_count = []

	for i in range(length):
		l_score[i] -= 1
		left_count = 0
		min_left = min(l_score)
		for j in range(length):
			if l_score[j] == min_left:
				left_count += l_count[j]

		r_score[i] -= 1
		right_count = 0
		min_right = min(r_score)
		for j in range(length):
			if r_score[j] == min_right:
				right_count += r_count[j]

		temp_score.append(min_left + min_right + 2)
		temp_count.append(left_count * right_count)
		l_score[i] += 1
		r_score[i] += 1

	score[node] = temp_score
	solution_count[node] = temp_count

	print('Before score :', l_score, r_score)
	print('After score :', temp_score)
	print('Before count :', l_count, r_count)
	print('After count :', temp_count)

def initialize_internal_nodes(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'postorder'):
		initialize_score_count(nonterminal)



def main():
	if len(sys.argv) == 3:
		INPUT_TREE_FILE = os.path.abspath(sys.argv[1])
		OUTPUT_FILE = os.path.abspath(sys.argv[2])
		# raise IndexError("Usage: python3 tnet.py [input phylogeny file] [desired output file]")
	else:
		INPUT_TREE_FILE = 'test/RAxML_rootedTree.4'
		OUTPUT_FILE = 'output.tnet'

	input_tree = initialize_tree(INPUT_TREE_FILE)
	initialize_leaf_nodes(input_tree)
	initialize_internal_nodes(input_tree)



if __name__ == "__main__": main()
