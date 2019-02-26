from __future__ import division
import sys, glob
sys.path.append('/scicore/home/neher/GROUP/data/Enterovirus_D68_analysis/SVVC/src')
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import defaultdict
#import seaborn as sns

from create_allele_counts import load_allele_counts
from coverage_consensus_diversity import coverage, consensus
from minor_variant import trim_ac

def running_average(a, ws=10):
    return np.convolve(a, np.ones(ws, dtype=float)/ws, mode='same')


plt.ion()
if __name__ == '__main__':
    sample_location = 'samples_by_strain/one_sample_per_patient/'

    freqs = {}
    diversity = {}
    cov = {}
    major_freq = {}
    samples = glob.glob(sample_location+'*')
    ref='KX675261.1'
    cutoff=0.01
    fs=16
    L=7333
    cross_sectional = np.zeros((6,L), dtype=float)
    coinf = ["EVD68_SWE_045_160831_NFLG", "EVD68_SWE_046_160904_NFLG", "EVD68_SWE_007_140908_NFLG"]

    for sample in samples:
        ac,ins = load_allele_counts(sample)
        sname = sample.rstrip('/').split('/')[-1]
        cov[sname] = coverage(ac[0][1] )
        freqs[sname] = trim_ac(ac, n_states=5)[ref]
        diversity[sname] = 1 - np.sum(freqs[sname]**2, axis=0)
        major_freq[sname] = np.max(freqs[sname], axis=0)
        cross_sectional[np.argmax(ac[0][1], axis=0), np.arange(L)] += (cov[sname]>100)

    cross_sectional /= (1e-10 + cross_sectional.sum(axis=0))
    entropy = -np.sum(cross_sectional*np.log(cross_sectional+1e-10), axis=0)

    min_cov = 1000
    snames = sorted(list(diversity.keys()))
    diversity_smooth = {}
    for sname in snames:
        good_ind = (cov[sname]>min_cov)&(major_freq[sname]<1-cutoff)
        all_pos = np.arange(good_ind.shape[0])
        tmp = []
        for i in range(3):
            rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==i)
            tmp.append(running_average(diversity[sname][rf]*good_ind[rf], ws=30))
        diversity_smooth[sname] = np.array(tmp)

    for i in range(3):
        plt.figure()
        for sname in diversity_smooth:
            if sname not in coinf:
                plt.plot(diversity_smooth[sname][i], c='C%d'%i)
        plt.yscale('log')
        plt.ylim([0.0001, 0.03])

    average_diversity = np.ma.array([np.ma.array((1-major_freq[s])*(major_freq[s]<1-cutoff), mask=(cov[s]<min_cov))
                                     for s,x in diversity.items() if s not in coinf]).max(axis=0)

    proteins = {'VP4':[699,906], 'VP2':[906,1650], 'VP3':[1650,2355], 'VP1':[2355,3282],
                '2A':[3282,3723], '2B':[3723, 4020], '2C':[4020,5010], '3A':[5010,5277],
                #'3B':[5277,5343],
                '3C':[5343, 5892], '3D':[5892,7263]}
    loops = {'BC':[89*3,103*3], 'DE':[140*3, 153*3]}
    for p,pos in proteins.items():
        proteins[p]=(np.array(pos)-699)/3
    for p,pos in loops.items():
        loops[p]=(np.array(pos))/3

    fig, axs = plt.subplots(3,2, sharex='col', figsize=(12,6))
    cp_labels = {0:'1st', 1:'2nd', 2:'3rd'}
    for ci, annos in [[0, proteins], [1,loops]]:
        for i in range(3):
            if ci==0:
                rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==i)
            else:
                rf = (all_pos>2355)&(all_pos<3282)&((all_pos-699)%3==i)

            #axs[0].plot(running_average(entropy[rf],ws=30), label='pos '+str(i+1))
            #axs[i].plot(entropy[rf],'o', label='pos '+str(i+1),)
            #axs[i].plot(running_average(average_diversity[rf],ws=30), label='pos '+str(i+1))
            axs[i, ci].plot(average_diversity[rf], '-o', c='C%d'%(i+1),
                            label=cp_labels[i]+' position')
            axs[i, ci].set_yscale('log')
            axs[i, ci].set_ylim([0.006, 3])
            axs[i, ci].text(1000,0.3 if i==1 else 1, cp_labels[i]+' codon position', fontsize=fs)
            axs[i, ci].tick_params(labelsize=fs*0.8)
            if i==1:
                if ci==0:
                    axs[i,ci].set_ylabel('iSNV frequency', fontsize=fs)
                pi = 0
                y1,y2= (0.003, 2) if ci else (1,2)
                for p, pos in sorted(annos.items(), key=lambda x:x[1][0]):
                    r = Rectangle((pos[0], y1),
                                  pos[1]-pos[0],
                                  y2-y1,
                                  facecolor=[0.6+0.2*(pi%2)] * 3,
                                  edgecolor='k',
                                  label=p)
                    pi+=1
                    xt = pos[0] + 0.5 * (pos[1]-pos[0])
                    yt = y2*0.53

                    axs[i,ci].add_patch(r)
                    axs[i,ci].text(xt, yt,
                        p,
                        color='k',
                        fontsize=fs*0.7,
                        ha='center')



    axs[-1,0].set_xlabel('polyprotein', fontsize=fs)
    axs[-1,1].set_xlabel('VP1', fontsize=fs)
    plt.tight_layout()
    for fmt in ['.pdf', '.png']:
        plt.savefig('figs/within_host_diversity_Fig4'+fmt)
        
