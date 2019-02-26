from __future__ import division
import sys, glob
sys.path.append('/home/richard/scicore/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
sys.path.append('/scicore/home/neher/GROUP/data/2017_Karolinska_EV-D68/SVVC/src')
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
#import seaborn as sns

from create_allele_counts import load_allele_counts
from coverage_consensus_diversity import coverage, consensus
from minor_variant import trim_ac

def peptide_to_nuc(pos):
    return 699 + (pos-1)*3

def load_ctl_prediction(ctl_fname):
    import pandas as pd
    pred = pd.read_csv(ctl_fname)
    regions = []
    for ri, row in pred.iterrows():
        nuc1,nuc2 = peptide_to_nuc(row.start),peptide_to_nuc(row.end+1)
        regions.append([row.allele, nuc1, nuc2, row.peptide, row.percentile_rank])

    return regions

def div_func(a):
    return a.mean()
#    return (a>0.01).mean()

plt.ion()
if __name__ == '__main__':

    freqs = {}
    diversity = {}
    cov = {}
    samples = glob.glob('mapped_data/1*')
    ctl_predictions = glob.glob('CTL_predictions/1*csv')
    predictions_by_sample = defaultdict(list)
    for ctl_fname in ctl_predictions:
        sname = ctl_fname.rstrip('/').split('/')[-1].split('_')[0]
        predictions_by_sample[sname].append(ctl_fname)

    ref='KX675261.1'
    comparison = [[], [],[]]
    for sname in predictions_by_sample:
        ctls = []
        for ctl_fname in predictions_by_sample[sname]:
            ctls.extend(load_ctl_prediction(ctl_fname))
        ctls.sort(key=lambda x:x[-1])

        ac,ins = load_allele_counts('mapped_data/'+sname)

        cov[sname] = coverage(ac[0][1] )
        freqs[sname] = trim_ac(ac, n_states=5)
        diversity[sname] = {ref:1 - np.sum(x**2, axis=0) for ref, x in freqs[sname].items()}
        CTL_mask = np.zeros_like(diversity[sname][ref], dtype=bool)
        good_ind = cov[sname]>1000
        all_pos = np.arange(good_ind.shape[0])
        for cutoff in [50]:
            for region in ctls[:cutoff]:
                CTL_mask[region[1]:region[2]]=True

            for i in range(3):
                rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==i)
                epi_div = div_func(diversity[sname][ref][CTL_mask&good_ind&rf])
                nonepi_div = div_func(diversity[sname][ref][(~CTL_mask)&good_ind&rf])
                print(epi_div, nonepi_div)
                if np.isnan(epi_div):
                    continue
                comparison[i].append((epi_div, nonepi_div))
    comparison = np.array(comparison)
    dcomp = comparison[:,:,1]-comparison[:,:,0]
    print(np.median(dcomp, axis=1))
    print(comparison.mean(axis=1))
