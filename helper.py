# Helper and other package testing for TNet in python3
# Saurav Dhar

from Bio import Phylo
import random
# from Bio.Phylo import NewickIO

easy_colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']


def calculate_score(left, right):
	result = []
	for i in range(len(left)):
        # Find min-cost change from node on the left
		temp_left = left.copy()
		temp_left[i] -= 1
		min_left = min(temp_left)

		# Find min-cost change from node on the right
		temp_right = right.copy()
		temp_right[i] -= 1
		min_right = min(temp_right)

		result.append(min_left + min_right + 2)

	return result


def calculate_child_state(current, left, right):
	result = []
	index = current.index(min(current))
	temp_left = left.copy()
	temp_left[index] -= 1
	min_left = min(temp_left)
	result.append(hosts[temp_left.index(min_left)])

	temp_right = right.copy()
	temp_right[index] -= 1
	min_right = min(temp_right)
	result.append(hosts[temp_right.index(min_right)])

	return result


rooted_tree = Phylo.read('test/RAxML_rootedTree.2', 'newick')

rooted_tree.rooted = True
print(rooted_tree.is_bifurcating())
# print(rooted_tree)


score = {}
child_state = {}
hosts = []
colors = []
all_clades = []
transmission_edges = []
r = lambda: random.randint(0,255)

print('Terminals: ', len(rooted_tree.get_terminals()))
print('Nonterminals: ', len(rooted_tree.get_nonterminals()))

for terminal in rooted_tree.get_terminals():
	terminal.name = terminal.name.split('_')[0]
	hosts.append(terminal.name)
	# print(terminal.name)

# Phylo.draw_ascii(rooted_tree)

# print(len(hosts))
hosts = list(set(hosts))
print(hosts)
print('Total hosts: ', len(hosts))

for host in hosts:
	color = '#%02X%02X%02X' % (r(),r(),r())
	colors.append(color)

# print('Colors: ', colors)

for terminal in rooted_tree.get_terminals():
	temp = []
	for host in hosts:
		if host==terminal.name:
			temp.append(0)
		else:
			temp.append(9999999999)

	score[terminal] = temp
	# print(terminal)

all_clades_size = len(rooted_tree.get_terminals()) + len(rooted_tree.get_nonterminals())

for nonterminal in rooted_tree.get_nonterminals(order = 'postorder'):
	# print('Calculating :',score[nonterminal.clades[0]], score[nonterminal.clades[1]])
	score[nonterminal] = calculate_score(score[nonterminal.clades[0]], score[nonterminal.clades[1]])
	child_state[nonterminal] = calculate_child_state(score[nonterminal], score[nonterminal.clades[0]], score[nonterminal.clades[1]])
	print(child_state[nonterminal])
	# print('After calcu :',score[nonterminal.clades[0]], score[nonterminal.clades[1]])

# Set the name of internal clades from score
rooted_tree.root.name = hosts[score[rooted_tree.root].index(min(score[rooted_tree.root]))]
rooted_tree.root.color = colors[hosts.index(rooted_tree.root.name)]
for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
	if not nonterminal.clades[0].is_terminal():
		nonterminal.clades[0].name = child_state[nonterminal][0]
		nonterminal.clades[0].color = colors[hosts.index(nonterminal.clades[0].name)]

	if not nonterminal.clades[1].is_terminal():
		nonterminal.clades[1].name = child_state[nonterminal][1]
		nonterminal.clades[1].color = colors[hosts.index(nonterminal.clades[1].name)]
	


# Get the transmission edges from the lebeled tree
for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
	if nonterminal.name != nonterminal.clades[0].name:
		print(score[nonterminal], score[nonterminal.clades[0]])
		print(nonterminal.name, nonterminal.clades[0].name)
		transmission_edges.append([nonterminal.name, nonterminal.clades[0].name])
	if nonterminal.name != nonterminal.clades[1].name:
		transmission_edges.append([nonterminal.name, nonterminal.clades[1].name])
		print(score[nonterminal], score[nonterminal.clades[1]])
		print(nonterminal.name, nonterminal.clades[1].name)


print('Transmission count:', len(transmission_edges), transmission_edges)
print('All clades: ', all_clades_size)
# print(all_clades)


print('Root: ', score[rooted_tree.root])
print('The minimum parsimony cost is:', min(score[rooted_tree.root]), 'with root:', rooted_tree.root.name)
# rooted_tree.ladderize()
# Phylo.draw(rooted_tree)


