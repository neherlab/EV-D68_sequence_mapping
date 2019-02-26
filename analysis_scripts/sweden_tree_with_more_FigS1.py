import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from Bio import Phylo

#old code:
# iqtree -s all_sequences_plus_dual_aln.fasta -ninit 2 -n 2 -me 0.05 -m GTR

# Trees are created by 'make_alignments.py' in Enterovirus_D68_analysis folder
# This uses the tree from the alignment that includes dually infected samples
treefile = 'alignments/subset-all-plus-swed-double_Ned.nwk'

T = Phylo.read(treefile, 'newick')

#T.root_at_midpoint()
T.root_with_outgroup("CAN12-106","EV-D68/Homosapiens/USA/C6840/2009")
T.ladderize()

#labels = {x.name:'_'.join([x.name.split('_')[i] for i in [1,2,3,5] if len(x.name.split('_'))>i]) for x in T.get_terminals()}
labels = {x.name: x.name.replace('EVD68_','').replace('_NFLG','') for x in T.get_terminals()}
labels[None]=''

#matplotlib.rc('font', size=7)

ib = T.common_ancestor("CAN12-106","2013-1017-26")
ib.name = "A1"
labels[ib.name] = ib.name

ib = T.common_ancestor("T106/FtJacksonSouthCarolinaUSA/2002", "EV-D68/Homosapiens/USA/C6840/2009")
ib.name = "A2"
labels[ib.name] = ib.name

ib = T.common_ancestor("4310902042_RD2", "CAN14-235")
ib.name = "B2"
labels[ib.name] = ib.name

ib = T.common_ancestor("EV-D68/Homosapiens/USA/U4489/2012", "US/CO/14-93")
ib.name = "B1"
labels[ib.name] = ib.name

ib = T.common_ancestor("TW-00909-2014", "EVD68_SWE_035_160829_NFLG")
ib.name = "B3"
labels[ib.name] = ib.name

def label_colors(x):
	print(x)
	if '007_14' in x:
		return 'C1'
	elif '045_16' in x:
		return 'C2'
	elif '046_16' in x:
		return 'C3'
	elif 'SWE_' in x:
		return 'C0'
	else:
		return 'k'

fig,ax = plt.subplots(1,1, figsize=(12,12))
Phylo.draw(T, label_func=lambda x:labels[x.name], axes=ax, label_colors=label_colors, show_confidence=False)
ax.set_axis_off()
plt.tight_layout()

plt.savefig('figs/sweden_and_more_tree_FigS1.pdf')
