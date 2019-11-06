#!/usr/bin/python3
# Here is the implementstion of TNet in python3
# Saurav Dhar

# Library imports
from Bio import Phylo
import numpy as np
import sys, os

# Global variables
score = {}
left_score = {}
right_score = {}
solution_count = {}
hosts = []
transmission_edges = []




def initialize_tree(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True
	# print('Terminals:', len(input_tree.get_terminals()))
	# print('Nonterminals:', len(input_tree.get_nonterminals()))
	if not input_tree.is_bifurcating():
		raise IndexError("Input tree is not bifurcating.")
	# Phylo.draw_ascii(rooted_tree)
	return input_tree

def initialize_leaf_nodes(rooted_tree):
	temp_host = []
	for terminal in rooted_tree.get_terminals():
		terminal.name = terminal.name.split('_')[0]
		temp_host.append(terminal.name)

	global hosts
	hosts = list(set(temp_host))
	hosts.sort()
	# print('Total hosts: ', len(hosts), hosts)

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
	temp_left = []
	temp_right = []
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
		temp_left.append(min_left)
		temp_right.append(min_right)
		temp_count.append(left_count * right_count)
		l_score[i] += 1
		r_score[i] += 1

	score[node] = temp_score
	left_score[node] = temp_left
	right_score[node] = temp_right
	solution_count[node] = temp_count

	# print('Before score :', l_score, r_score)
	# print('Left:', temp_left, 'Right:', temp_right)
	# print('After score :', temp_score)
	# print('Before count :', l_count, r_count)
	# print('After count :', temp_count)
	# print('=========================')

def initialize_internal_nodes(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'postorder'):
		initialize_score_count(nonterminal)

def get_host_from_count(count):
	probs = [float(i)/sum(count) for i in count]
	ch = np.random.choice(len(probs), p=probs)
	return hosts[ch]

def choose_root_host(root_node):
	probs = []
	min_score = min(score[root_node])
	for i in range(len(score[root_node])):
		if score[root_node][i]==min_score:
			probs.append(solution_count[root_node][i])
		else:
			probs.append(0)

	# print('Root', probs)
	# print('Root score', score[root_node])
	return get_host_from_count(probs)

def choose_internal_node_host(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		# print(score[nonterminal])
		# print(solution_count[nonterminal])
		index = hosts.index(nonterminal.name)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			l_count = solution_count[nonterminal.clades[0]].copy()
			l_score[index] -= 1
			for i in range(len(l_score)):
				if l_score[i] != left_score[nonterminal][index]:
					l_count[i] = 0

			nonterminal.clades[0].name = get_host_from_count(l_count)

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			r_count = solution_count[nonterminal.clades[1]].copy()
			r_score[index] -= 1
			for i in range(len(r_score)):
				if r_score[i] != right_score[nonterminal][index]:
					r_count[i] = 0

			nonterminal.clades[1].name = get_host_from_count(r_count)


def get_transmission_edges(rooted_tree):
	edges = []
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		if nonterminal.name != nonterminal.clades[0].name:
			edges.append([nonterminal.name, nonterminal.clades[0].name])
		if nonterminal.name != nonterminal.clades[1].name:
			edges.append([nonterminal.name, nonterminal.clades[1].name])

	return edges

def write_transmission_edges(file, source, edges):
	result = open(file, 'w+')
	result.write('{}\t{}\n'.format('None', source))
	for edge in edges:
		result.write('{}\t{}\n'.format(edge[0], edge[1]))

	result.close()

def main():
	if len(sys.argv) == 3:
		INPUT_TREE_FILE = os.path.abspath(sys.argv[1])
		OUTPUT_FILE = os.path.abspath(sys.argv[2])
	else:
		raise IndexError("Usage: python3 tnet.py [input phylogeny file] [desired output file]")

	input_tree = initialize_tree(INPUT_TREE_FILE)
	initialize_leaf_nodes(input_tree)
	initialize_internal_nodes(input_tree)
	input_tree.root.name = choose_root_host(input_tree.root)
	choose_internal_node_host(input_tree)
	transmission_edges = get_transmission_edges(input_tree)
	write_transmission_edges(OUTPUT_FILE, input_tree.root.name, transmission_edges)

	# print('Transmission count:', len(transmission_edges), transmission_edges)
	print('The minimum parsimony cost is:', min(score[input_tree.root]), 'with root:', input_tree.root.name)

if __name__ == "__main__": main()
