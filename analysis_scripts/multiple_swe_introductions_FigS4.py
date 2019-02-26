import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from Bio import Phylo


# IQ-tree command by Richard used to generate tree for 2014-2016 paper.
'''
Host:    ran13 (AVX2, FMA3, 15 GB RAM)
Command: iqtree -s aligned_genome.fasta -bb 1000 -m HKY
Seed:    935247 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Tue Feb  5 15:50:43 2019
Kernel:  AVX+FMA - 1 threads (4 CPU cores detected)
'''
#This tree was:
# /home/richard/Projects/nextstrain/enterovirus/bootstrap_tree/aligned_genome.fasta.treefile

#For future runs, upload the alignment from the Nextstrain run to 'alignments'
#Run iqtree on aligned_genome.fasta to get a tree:
#iqtree -s alignments/aligned_genome.fasta -bb 1000 -m HKY

treefile = 'alignments/aligned_genome.fasta.treefile' 

sweden2016_color = "#33AA33"
other_color = "#666666"
cutoff = 0.9

T = Phylo.read(treefile, 'newick')
T.ladderize()
Phylo.draw(T)


nodes_to_collapse = []
for n in T.find_clades():
	if n.confidence and n.confidence*0.01<cutoff:
		nodes_to_collapse.append(n)

for n in nodes_to_collapse:
	T.collapse(n)

#These nodes should be checked in future. For example, now the second one would be: EVD68/Homosapiens/VNM/VN11/2015
tmp_nodes = [n for n in T.get_terminals() if n.name in ['EVD68_SWE_031_160822_NFLG', 'EVD68_Homosapiens_VNM_VN11_2015']]
sweden2016 = T.common_ancestor(tmp_nodes)


for n in T.find_clades(order='postorder'):
	if n.is_terminal():
		n.mono_sweden = ('SWE_' in n.name) and ('_16' in n.name)
	else:
		n.mono_sweden = all([c.mono_sweden for c in n])

for n in T.find_clades(order='preorder'):
	n.color = sweden2016_color if n.mono_sweden else other_color
	print(n.mono_sweden, n.color)

for n in T.get_nonterminals():
	for c in n:
		if c.mono_sweden and (not n.mono_sweden):
			print('putative introduction', c.name, len(c.get_terminals()))

fig = plt.figure(figsize=(10,13))
ax = plt.subplot(111)

Phylo.draw(sweden2016, label_colors=lambda x:sweden2016_color if 'SWE' in x else other_color, axes=ax)
plt.axis('off')
plt.tight_layout()
plt.savefig('figs/multiple_swe_introductions_FigS4.pdf')
