#!/usr/bin/python3
"""
Copyright (C) 2019 Saurav Dhar (saurav.dhar@uconn.edu), Chengchen Zhang,
Ion Mandoiu (ion.mandoiu@uconn.edu), and Mukul S. Bansal
(mukul.bansal@uconn.edu).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Library imports
from Bio import Phylo
import numpy as np
import sys, os
import argparse

# Global variables
score = {}
left_score = {}
right_score = {}
solution_count = {}
hosts = []
transmission_edges = []
rand_seed = None
flag_max_prob = None
flag_equal_prob = None


def initialize_tree(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True

	if not input_tree.is_bifurcating():
		raise IndexError("Input tree is not bifurcating.")

	return input_tree

def initialize_leaf_nodes(rooted_tree):
	temp_host = []
	for terminal in rooted_tree.get_terminals():
		terminal.name = terminal.name.split('_')[0]
		if terminal.name not in temp_host:
			temp_host.append(terminal.name)

	global hosts
	hosts = temp_host.copy()

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

def initialize_internal_nodes(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'postorder'):
		initialize_score_count(nonterminal)

def get_host_from_count(count):
	if flag_equal_prob:
		for i in range(len(count)):
			if count[i] != 0:
				count[i] = 1
	elif flag_max_prob:
		max_count = max(count)
		for i in range(len(count)):
			if count[i] != max_count:
				count[i] = 0

	probs = [float(i)/sum(count) for i in count]
	np.random.seed(rand_seed)
	ch = np.random.choice(len(probs), p=probs)
	return hosts[ch]

def choose_root_host(root_node):
	probs = []
	min_score = min(score[root_node])
	for i in range(len(score[root_node])):
		if score[root_node][i] == min_score:
			probs.append(solution_count[root_node][i])
		else:
			probs.append(0)

	return get_host_from_count(probs)

def choose_internal_node_host(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
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

def choose_internal_node_host_with_bias(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		index = hosts.index(nonterminal.name)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			l_count = solution_count[nonterminal.clades[0]].copy()
			if min(l_score) == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			elif min(l_score) + 1 == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			else:
				countTotal = 0
				for i in range(len(l_score)):
					if l_score[i] == min(l_score):
						countTotal += l_count[i]

				r = np.random.randint(countTotal)
				for i in range(len(l_score)):
					if l_score[i] == min(l_score):
						r -= l_count[i]
						if r <= 0:
							nonterminal.clades[0].name = hosts[i]
							break

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			r_count = solution_count[nonterminal.clades[1]].copy()
			if min(r_score) == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			elif min(r_score) + 1 == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			else:
				countTotal = 0
				for i in range(len(r_score)):
					if r_score[i] == min(r_score):
						countTotal += r_count[i]

				r = np.random.randint(countTotal)
				for i in range(len(r_score)):
					if r_score[i] == min(r_score):
						r -= r_count[i]
						if r <= 0:
							nonterminal.clades[1].name = hosts[i]
							break

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

def write_info_file(file, rooted_tree):
	result = open(file + '.info', 'w+')
	result.write('Hosts list: {}\n'.format(hosts))
	result.write('Root ID: {}\n'.format(rooted_tree.root.name))
	result.write('The minimum parsimony cost: {}\n'.format(min(score[rooted_tree.root])))
	result.write('Root score list: {}\n'.format(score[rooted_tree.root]))
	result.write('Root assignment count list: {}\n'.format(solution_count[rooted_tree.root]))

	result.write('\nChanges:\n')
	edges = get_transmission_edges(rooted_tree)
	for edge in edges:
		result.write('{} -> {}\n'.format(edge[0], edge[1]))

	result.write('\nNewick Format Tree:\n')
	Phylo.write([rooted_tree], result, 'newick')
	result.close()

def read_parser_args(args):
	global rand_seed
	rand_seed = args.seed
	global flag_max_prob
	flag_max_prob = args.maxprob
	global flag_equal_prob
	flag_equal_prob = args.equalprob

def main():
	parser = argparse.ArgumentParser(description='Process TNet arguments.')
	parser.add_argument('INPUT_TREE_FILE', action='store', type=str, help='input file name')
	parser.add_argument('OUTPUT_FILE', action='store', type=str, help='output file name')
	parser.add_argument('-sd', '--seed', default=None, type=int, help='random number generator seed')
	parser.add_argument('-rs', '--randomsampling', default=False, action="store_true", help='choose internal node host with probability proportional to solutions')
	parser.add_argument('-mx', '--maxprob', default=False, action="store_true", help='choose internal node host with max probability')
	parser.add_argument('-eq', '--equalprob', default=False, action="store_true", help='choose internal node host with all equal probability')
	parser.add_argument('-info', '--info', default=False, action="store_true", help='write info file')
	parser.add_argument('--version', action='version', version='%(prog)s 1.1')
	args = parser.parse_args()

	read_parser_args(args)
	input_tree = initialize_tree(args.INPUT_TREE_FILE)
	initialize_leaf_nodes(input_tree)
	initialize_internal_nodes(input_tree)
	input_tree.root.name = choose_root_host(input_tree.root)
	if args.randomsampling:
		choose_internal_node_host(input_tree)
	else:
		choose_internal_node_host_with_bias(input_tree)

	transmission_edges = get_transmission_edges(input_tree)
	write_transmission_edges(args.OUTPUT_FILE, input_tree.root.name, transmission_edges)

	if args.info:
		write_info_file(args.OUTPUT_FILE, input_tree)

	print('The minimum parsimony cost is:', min(score[input_tree.root]), 'with root:', input_tree.root.name)

if __name__ == "__main__": main()
