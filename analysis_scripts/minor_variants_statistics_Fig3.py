from __future__ import division
import sys, glob
sys.path.append('/scicore/home/neher/GROUP/data/Enterovirus_D68_analysis/SVVC/src')
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#import seaborn as sns

from create_allele_counts import load_allele_counts
from coverage_consensus_diversity import coverage, consensus
from minor_variant import trim_ac
from helpers import add_panel_label

plt.ion()
if __name__ == '__main__':
    sample_location = 'samples_by_strain/one_sample_per_patient/'
    #sample_location = 'samples_by_strain/EVD68_SWE_0'

    freqs = {}
    major_freqs = {}
    cov = {}
    samples = glob.glob(sample_location+'*')

    for sample in samples:
        ac,ins = load_allele_counts(sample)
        sname = sample.rstrip('/').split('/')[-1]
        cov[sname] = coverage(ac[0][1] )
        freqs[sname] = trim_ac(ac, n_states=5)
        major_freqs[sname] = {ref:np.max(x, axis=0) for ref, x in freqs[sname].items()}

    #cutoffs = [0.01, 0.03, 0.1]
    cutoffs = [0.0005,0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1]
    display_cutoffs = [0.01, 0.03, 0.1]
    min_cov = 1000
    n_minor = []
    snames = sorted(list(major_freqs.keys()))
    ref='KX675261.1'
    for sname in snames:
        good_ind = cov[sname]>min_cov
        tmp = []
        tmp.append([np.sum((major_freqs[sname][ref]<1.0-c)&good_ind) for c in cutoffs])
        all_pos = np.arange(good_ind.shape[0])
        for i in range(3):
            rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==i)
            tmp.append([np.sum((major_freqs[sname][ref]<1.0-c)&(good_ind&rf)) for c in cutoffs])
        n_minor.append(tmp)

    n_minor = np.array(n_minor)

    fs=16
    bw = 2
    lst = '-', '--', ':', '-.'

    sorted_ii = np.argsort(n_minor[:,-1,-3])
    nii=5
    print(np.array(snames)[sorted_ii[-nii:]])
    print(n_minor[sorted_ii[-nii:]])

    cfrac = []
    for ci, c in enumerate(cutoffs):
        tmp = []
        good_ind = n_minor[:,1:,ci].sum(axis=1)>0
        for ii in [1,2,3]:
            dnds = (1.0*n_minor[:,ii,ci]/n_minor[:,1:,ci].sum(axis=1))[good_ind]
            tmp.append(dnds.mean())
        print(c,tmp, good_ind.sum())
        cfrac.append(tmp)
    cfrac = np.array(cfrac)


    fs=24
    fig, axs = plt.subplots(1,2, figsize=(12,6))
    ci = cutoffs.index(0.03)
    cp_labels = {1:'1st', 2:'2nd', 3:'3rd'}
    for i in range(4):
        axs[1].plot(sorted(n_minor[:,i,ci]), n_minor.shape[0] - np.arange(n_minor.shape[0]),
                 '-' if i else '-o', label=cp_labels[i] if i else "all",
                 lw=3, c='C%d'%i)
                 #label="frequency cutoff "+str(c), lw=3, ls=lst[i], c='C%d'%ci)

    axs[1].set_yscale('log')
    axs[1].set_xscale('log')
    axs[1].set_ylabel('number of samples', fontsize=fs)
    axs[1].set_xlabel('number of minor variants', fontsize=fs)
    axs[1].tick_params(labelsize=0.8*fs)
    axs[1].legend(fontsize=fs, loc=1, frameon=True, edgecolor='k')
    add_panel_label(axs[1], 'B', fs, x_offset=-0.18)

    for i in range(3):
        axs[0].plot(cutoffs, cfrac[:,i], label='codon pos %d'%(i+1), c='C%d'%(i+1), lw=3)

    axs[0].set_ylabel('fraction of variable sites', fontsize=fs)
    axs[0].set_xlabel('frequency cutoff', fontsize=fs)
    axs[0].set_xscale('log')
    axs[0].set_ylim(0.0, 0.7)
    axs[0].tick_params(labelsize=0.8*fs)
    add_panel_label(axs[0], 'A', fs, x_offset=-0.18)
    plt.tight_layout()
    for fmt in ['.pdf', '.png']:
        plt.savefig('figs/minor_variants_Fig3'+fmt)
    


