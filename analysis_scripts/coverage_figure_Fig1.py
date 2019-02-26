from __future__ import division
import sys, glob
#sys.path.append('/home/richard/scicore/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
sys.path.append('/scicore/home/neher/GROUP/data/Enterovirus_D68_analysis/SVVC/src')
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
#import seaborn as sns

from create_allele_counts import load_allele_counts, nuc_alpha
from coverage_consensus_diversity import coverage, consensus, get_fragment_boundaries
from minor_variant import trim_ac
from helpers import name_translations


if __name__ == '__main__':
    #ntt = name_translations('name_translation_table.tsv')
    #labid = '16CA403716'
    #name = ntt[labid]
    sample_location = 'samples_by_strain/'
    name = "EVD68_SWE_029_160904_NFLG"
    fs=16

    ac, ins = load_allele_counts(sample_location+name+"/")
    primer_boundaries = get_fragment_boundaries('primers.csv', ac)
    cov = coverage(ac[0][1])
    ref='KX675261.1'

    plt.figure(figsize=(8,4))
    plt.plot(cov, lw=3)
    for p in primer_boundaries[ref]:
        y = 50 if int(p[1])%2 else 70
        plt.plot([primer_boundaries[ref][p]['start'], primer_boundaries[ref][p]['end']],[y,y], lw=7, c=(0.7, 0.7, 0.7))

    plt.xlabel('position in genome', fontsize=fs)
    plt.ylabel('coverage', fontsize=fs)
    plt.ylim(30,10000)
    plt.yscale('log')
    plt.tick_params(labelsize=0.8*fs)
    plt.title(name, fontsize=fs)
    plt.tight_layout()
    for fmt in ['pdf', 'png']:
        plt.savefig('figs/coverage_example_Fig1.'+fmt)
