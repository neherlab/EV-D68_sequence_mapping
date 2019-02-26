import sys
#sys.path.append('/home/richard/scicore/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
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
    #Sample A is SWE_027
    #Sample B is SWE_037
    nameKey = {'A':'SWE_027', 'B':'SWE_037'}

    freqs = {}
    major_freqs = {}
    cov = {}
    for sample in ['QC_sample-A_run1', 'QC_sample-A_run2', 'QC_sample-B_run1', 'QC_sample-B_run2']:    
        ac,ins = load_allele_counts(sample_location + sample + "/")
        cov[sample] = coverage(ac[0][1] )
        freqs[sample] = trim_ac(ac, n_states=5)
        major_freqs[sample] = {ref:np.max(x, axis=0) for ref, x in freqs[sample].items()}


    fs=24
    fig, axs = plt.subplots(1,2, figsize=(12,6))
    axs[0].plot([0.001, 1.0], [0.001, 1.0], c='k')
    min_cov=2000
    for s in 'AB':
        for ref in major_freqs['QC_sample-'+s+'_run1']:
            good_ind = (cov['QC_sample-'+s+'_run1']>min_cov)&(cov['QC_sample-'+s+'_run2']>min_cov)
            axs[0].scatter(1.0-major_freqs['QC_sample-'+s+'_run1'][ref][good_ind],
                        1.0-major_freqs['QC_sample-'+s+'_run2'][ref][good_ind], label=nameKey[s])
                        
    axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    axs[0].set_ylim(0.001, .5)
    axs[0].set_xlim(0.001, .5)
    axs[0].set_ylabel('minor freqencies in run 2', fontsize=fs)
    axs[0].set_xlabel('minor freqencies in run 1', fontsize=fs)
    axs[0].tick_params(labelsize=0.8*fs)
    axs[0].legend(fontsize=fs, frameon=True, edgecolor='k', loc=2)
    add_panel_label(axs[0], 'A', fs, x_offset=-0.18)

    all_pos = np.arange(good_ind.shape[0])
    exclude = [690, 694]
    bins = np.logspace(-4,-2,11)
    bc=0.5*(bins[:-1]+bins[1:])
    lc = 0
    for s in 'AB':
        for ref in major_freqs['QC_sample-'+s+'_run1']:
            for ri, r in enumerate(['_run1', '_run2']):
                lc+=1
                good_ind = cov['QC_sample-'+s+r]>min_cov
                good_ind[exclude]=False
                counts = np.bincount(np.searchsorted(bins, 1.0-major_freqs['QC_sample-'+s+r][ref][good_ind]))[1:-1]
                axs[1].plot(bc[:len(counts)], 1.0*counts/counts.sum(),
                            label=nameKey[s]+", run %s"%(ri+1), lw=3, c="C%d"%lc)

                # !!! this commented out code has not been converted to new naming!!
                # there is essentially no between codon positions in the bulk
                # rf2 = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==2)
                # counts = np.bincount(np.searchsorted(bins, 1.0-major_freqs[r+'-'+s][ref][good_ind&rf2]))[1:-1]
                # axs[1].plot(bc[:len(counts)], 1.0*counts/counts.sum(),
                #             ls='--', lw=2, c="C%d"%lc)

    axs[1].set_xscale('log')
    axs[1].set_ylim(0.0, .49)
    axs[1].set_xlabel('minor freqency',fontsize=fs)
    axs[1].set_ylabel('fraction of sites',fontsize=fs)
    axs[1].tick_params(labelsize=0.8*fs)
    axs[1].legend(fontsize=fs*0.8, frameon=True, edgecolor='k')
    add_panel_label(axs[1], 'B', fs, x_offset=-0.18)

    plt.tight_layout()
    for fmt in ['.pdf', '.png']:
        plt.savefig('figs/QC_minor_variants_Fig2'+fmt)
