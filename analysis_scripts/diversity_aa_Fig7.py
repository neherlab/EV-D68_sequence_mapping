import numpy as np
import json
from Bio import Phylo
from collections import defaultdict
import matplotlib
matplotlib.use('agg')
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

with open('nextstrain_files/enterovirus_d68_genome_tree.json') as fh:
    T = json.load(fh)

current_state = defaultdict(lambda :'ancestral')
ancestral_states = {}
state_counts = defaultdict(lambda :defaultdict(int))
mutation_count = defaultdict(lambda: defaultdict(list))

pre_order_node_list = []
def pre_order(node):
    pre_order_node_list.append(node)
    if 'children'  in node and node["children"]:
        for c in node["children"]:
            c["parent"] = node
            pre_order(c)


post_order_node_list = []
def post_order(node):
    if 'children'  in node and node["children"]:
        for c in node["children"]:
            post_order(c)
    post_order_node_list.append(node)


post_order(T)
pre_order(T)
L=7500
fs=16
all_pos = np.arange(7500)
for n in pre_order_node_list[1:]:
    for prot in n["aa_muts"]:
        for m in n["aa_muts"][prot]:
            a,pos,d = m[0], int(m[1:-1]), m[-1]
            if d!='X':
                mutation_count[prot][pos].append(m)


proteins = {'VP4':[733,940], 'VP2':[940,1684], 'VP3':[1684,2389], 'VP1':[2389,3316],
            '2A':[3316,3757], '2B':[3757, 4054], '2C':[4054,5044], '3A':[5044,5311],
            '3B':[5311,5377],
            '3C':[5377, 5926], '3D':[5926,7297]}

loops = {'BC':[89*3,103*3], 'DE':[140*3, 153*3]}
for p,pos in proteins.items():
    proteins[p]=(np.array(pos)-730)/3
for p,pos in loops.items():
    loops[p]=(np.array(pos))/3

nmuts = {}
for p, pos in sorted(proteins.items(), key=lambda x:x[1][0]):
    nmuts[p] = np.array([len(mutation_count[p][i])
                    for i in range(int(pos[1]-pos[0]))])

fig, axs = plt.subplots(2,1, figsize=(8,6))
cp_labels = {0:'1st', 1:'2nd', 2:'3rd'}
pi = 0
ymax=0
ymin = -5
for ci, annos in [[0, proteins], [1,loops]]:
    for p, pos in sorted(annos.items(), key=lambda x:x[1][0]):
        if ci==0:
            print(p,pos, nmuts[p])
            axs[ci].plot(np.arange(pos[0], pos[1]), nmuts[p], '-', lw=2, c='C%d'%(pi%9+1))
        else:
            axs[ci].plot(nmuts['VP1'], '-', lw=2, c='C%d'%(pi%9+1))
        axs[ci].set_ylim([ymin,30])
        axs[ci].tick_params(labelsize=fs*0.8)
        axs[ci].set_ylabel('# of changes', fontsize=fs)
        y1,y2= (ymin, ymax) if ci else (ymin,ymax)
        r = Rectangle((pos[0], y1),
                      pos[1]-pos[0],
                      y2-y1,
                      facecolor='C%d'%(pi%9+1),
                      edgecolor='k',
                      alpha=0.3,
                      label=p)
        pi+=1
        xt = pos[0] + 0.5 * (pos[1]-pos[0])
        yt = y2-4

        axs[ci].add_patch(r)
        axs[ci].text(xt, yt,
            p,
            color='k',
            fontsize=fs*0.7,
            ha='center')



axs[0].set_xlabel('polyprotein', fontsize=fs)
axs[1].set_xlabel('VP1', fontsize=fs)
plt.tight_layout()
for fmt in ['.pdf', '.png']:
    plt.savefig('figs/diversity_aa_Fig7'+fmt)
