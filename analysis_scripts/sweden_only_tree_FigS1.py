import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from Bio import Phylo

#old code:
# iqtree -s all_sequences_plus_dual_aln.fasta -ninit 2 -n 2 -me 0.05 -m GTR

# Trees are created by 'make_alignments.py' in Enterovirus_D68_analysis folder
# This uses the tree from the alignment that includes dually infected samples
treefile = 'alignments/all_seq_inc_dual.nwk'

T = Phylo.read(treefile, 'newick')

T.root_at_midpoint()
T.ladderize()

labels = {x.name:'_'.join([x.name.split('_')[i] for i in [1,2,3,5] if len(x.name.split('_'))>i]) for x in T.get_terminals()}
labels[None]=''

def label_colors(x):
	print(x)
	if '007_14' in x:
		return 'C1'
	elif '045_16' in x:
		return 'C2'
	elif '046_16' in x:
		return 'C3'
	else:
		return 'k'


fig,ax = plt.subplots(1,1, figsize=(12,12))
Phylo.draw(T, label_func=lambda x:labels[x.name], axes=ax, label_colors=label_colors, show_confidence=False)
ax.set_axis_off()

plt.savefig('figs/sweden_only_tree_FigS1.pdf')
