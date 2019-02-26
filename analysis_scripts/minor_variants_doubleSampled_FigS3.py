import sys
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
    sample_location = 'samples_by_strain/'

    freqs = {}
    major_freqs = {}
    cov = {}
    samps = {'SWE_012': {'A':'EVD68_SWE_012_160911_NFLG', 'B':'EVD68_SWE_012_160913_NFLG'},
    'SWE_021': {'A':'EVD68_SWE_021_160828_NFLG', 'B':'EVD68_SWE_021_160904_NFLG'},
    'SWE_024': {'A':'EVD68_SWE_024_160830_NFLG', 'B':'EVD68_SWE_024_160831_NFLG'},
    'SWE_037': {'A':'EVD68_SWE_037_160829_NFLG', 'B':'EVD68_SWE_037_160830_NFLG'},
    'SWE_039': {'A':'EVD68_SWE_039_160831-1_NFLG', 'B':'EVD68_SWE_039_160831-2_NFLG'}}

    days = {'SWE_012': '2 days', 'SWE_021':'7 days', 'SWE_024':'1 day', 'SWE_037':'1 day', 'SWE_039':'0 days'}

    for pt in samps:
        for key in samps[pt]:
            ac,ins = load_allele_counts(sample_location+samps[pt][key])
            sample = pt+'-'+key
            cov[sample] = coverage(ac[0][1] )
            freqs[sample] = trim_ac(ac, n_states=5)
            major_freqs[sample] = {ref:np.max(x, axis=0) for ref, x in freqs[sample].items()}


##################################
##################################

    #positions not to plot
    #These positions are just before CDS, and show up in
    #4/5 of these samples....
    exclude = [690, 694]

    fs=12
    fig, axs = plt.subplots(2,3, figsize=(10,8))

    #make common axes...
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    plt.xlabel("minor frequencies in sample 1", fontsize=fs, labelpad=8)
    plt.ylabel("minor frequencies in sample 2", fontsize=fs, labelpad=8)

    for s in range(5):
        axs[s%2][s//2].plot([0.001, 1.0], [0.001, 1.0], c='k')
    min_cov=2000
    ns = 0
    for s in samps:
        for ref in major_freqs[s+'-A']:
            good_ind = (cov[s+'-A']>min_cov)&(cov[s+'-B']>min_cov)
            good_ind[exclude] = False
            axs[ns%2][ns//2].scatter(1.0-major_freqs[s+'-A'][ref][good_ind],
                        1.0-major_freqs[s+'-B'][ref][good_ind], label='patient '+s+'\ninterval: '+days[s],
                        c="C%d"%ns)
            axs[ns%2][ns//2].set_yscale('log')
            axs[ns%2][ns//2].set_xscale('log')
            axs[ns%2][ns//2].set_ylim(0.0005, .6)
            axs[ns%2][ns//2].set_xlim(0.0005, .6)
            axs[ns%2][ns//2].tick_params(labelsize=0.8*fs)
            axs[ns%2][ns//2].legend(fontsize=fs, frameon=True, edgecolor='k')
            ns+=1

    axs[1][2].set_axis_off()
    plt.tight_layout()
    for fmt in ['.pdf', '.png']:
        plt.savefig('figs/minorVar_doubleSampled_FigS3'+fmt)
