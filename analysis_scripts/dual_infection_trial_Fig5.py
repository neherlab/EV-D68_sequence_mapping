from __future__ import division
import sys, glob
sys.path.append('/scicore/home/neher/GROUP/data/Enterovirus_D68_analysis/SVVC/src')
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

from create_allele_counts import load_allele_counts, nuc_alpha
from coverage_consensus_diversity import coverage, consensus, get_fragment_boundaries
from minor_variant import trim_ac
from helpers import add_panel_label, name_translations



if __name__ == '__main__':
    #sample_location = 'samples_by_strain/'
    sample_location = "samples_by_strain2/"
    recalculate = True #False

    freqs = {}
    major_freqs = {}
    major_seqs = {}
    cov = {}
    #snames = ["EVD68_SWE_045_160831_NFLG", "EVD68_SWE_046_160904_NFLG", "EVD68_SWE_007_140908_NFLG"]
    #snames = ["EVD68_NLD_007_18XXXX_NFLG", "EVD68_BEL_009_18XXXX_NFLG"]
    snames = ["EVD68_BEL_009_18XXXX_NFLG"]    
    for sname in snames:
        ac,ins = load_allele_counts(sample_location+sname)
        primer_boundaries = get_fragment_boundaries('primers.csv', ac)
        cov[sname] = coverage(ac[0][1] )
        freqs[sname] = trim_ac(ac, n_states=5)
        major_freqs[sname] = {ref:np.max(x, axis=0) for ref, x in freqs[sname].items()}
        major_seqs[sname] = {ref:nuc_alpha[np.argmax(x, axis=0)] for ref, x in freqs[sname].items()}

    min_cov = 1000
    ref='KX675261.1'
    cp_labels = {0:'non coding', 1:'1st', 2:'2nd', 3:'3rd'}
    variable_sites = {}
    fig, axs = plt.subplots(len(snames)+1,2, figsize=(15,12),
                    gridspec_kw={"width_ratios":[2,1], "height_ratios":[1,0.1]}) # [1,1,1,0.1]})
    fs=16
    for sname, ax in zip(snames, axs[:,0]):
        print("Loading frequencies sample",sname)
        good_ind = cov[sname]>min_cov
        all_pos = np.arange(good_ind.shape[0])
        variable_sites[sname]=set()
        for cmin, cmax, a in [(0.03, 1.0, 0.9), (0.01,0.03, 0.6), (0.003,0.01, 0.3)]:
            for i in range(4):
                if i:
                    rf = (all_pos>699)&(all_pos<7266)&((all_pos-699)%3==(i-1))
                else:
                    rf = (all_pos<=699)|(all_pos>=7266)

                ind = good_ind&rf&(major_freqs[sname][ref]<1-cmin)&(major_freqs[sname][ref]>=1-cmax)
                if cmin>=0.01:
                    variable_sites[sname].update(np.where(ind)[0])

                ax.scatter(all_pos[ind], (1-major_freqs[sname][ref])[ind], c='C%d'%(i+1),
                            label=cp_labels[i] if cmax>0.9 else None, alpha=a)
                ax.set_yscale('log')
                ax.set_ylim([0.002, 1.0])

                #ax.text(100,0.5,ntt[sname][:-5], fontsize=fs)
                ax.text(100,0.5,sname[:-5], fontsize=fs)
                ax.tick_params(labelsize=0.8*fs)
    axs[0,0].legend(loc=1,ncol=2, fontsize=0.8*fs, frameon=True, edgecolor='k')
    axs[-2,0].set_xlabel('position in genome', fontsize=fs)
    axs[1,0].set_ylabel('frequency', fontsize=fs)

    for p in primer_boundaries[ref]:
        y = int(p[1])%2
        axs[-1,0].plot([primer_boundaries[ref][p]['start'], primer_boundaries[ref][p]['end']],[y,y], lw=7, c=(0.7, 0.7, 0.7))
    axs[-1,0].set_axis_off()
    axs[-1,0].set_ylim([-0.3,1.5])


    import gzip
    import cPickle as pickle
    from itertools import combinations
    if recalculate:
        for si,sname in enumerate(snames):
            print("Loading pair frequencies sample",sname)
            L = ac[0][1].shape[1]
            n = len(variable_sites[sname])
            pair_freq = np.zeros((n,n), dtype=float)
            freq = np.zeros(n, dtype=float)

            with gzip.open(sample_location+'%s/pair_counts.pkl.gz'%sname) as fh:
                ac, acc = pickle.load(fh)
            vsites = list(enumerate(sorted(variable_sites[sname])))
            pos_map = np.array([x[1] for x in vsites])
            for pi,pos in vsites:
                import ipdb; ipdb.set_trace()
                nuc = list(nuc_alpha).index(major_seqs[sname][ref][pos])
                tmp_cov = ac[0][1][:,pos].sum()
                if tmp_cov>100:
                    freq[pi] = 1.0*ac[0][1][nuc,pos]/tmp_cov
                else:
                    print(pos, tmp_cov)
                    freq[pi] = np.nan
                #freq[pi] = major_freqs[sname][ref][pos]

            pair_freq = np.zeros((n,n), dtype=float)
            p1_freq = np.zeros((n,n), dtype=float)
            p2_freq = np.zeros((n,n), dtype=float)
            print("Calculating LD",sname)
            for (pi1,p1),(pi2, p2)  in combinations(vsites, 2):
                posp = (p1,p2)
                if posp not in acc[0][1]:
                    continue

                nn = major_seqs[sname][ref][p1]+major_seqs[sname][ref][p2]
                if nn not in acc[0][1][posp]:
                    continue

                valid_gts = [x for x in acc[0][1][posp] if isinstance(x, str)]
                pair_cov = 1.0*sum([acc[0][1][posp][x] for x in valid_gts])
                if pair_cov>100:
                    pair_freq[pi1,pi2] = acc[0][1][posp][nn]/pair_cov
                    p1_freq[pi1,pi2] = sum([acc[0][1][posp][x] for x in valid_gts if x[0]==nn[0]])/pair_cov
                    p2_freq[pi1,pi2] = sum([acc[0][1][posp][x] for x in valid_gts if x[1]==nn[1]])/pair_cov

            pair_freq = pair_freq + pair_freq.T
            pair_freq = np.ma.array(pair_freq, mask=pair_freq==0)

            p1_freq = p1_freq + p1_freq.T
            p1_freq = np.ma.array(p1_freq, mask=p1_freq==0)

            p2_freq = p2_freq + p2_freq.T
            p2_freq = np.ma.array(p2_freq, mask=p2_freq==0)

            corr_matrix = (pair_freq - p1_freq*p2_freq)
            tmp = np.zeros_like(pair_freq)
            posLD = np.where(corr_matrix>0)
            negLD = np.where(corr_matrix<=0)
            tmp[posLD] = np.minimum(p1_freq*(1-p2_freq), (1-p1_freq)*p2_freq)[posLD]
            tmp[negLD] = np.minimum(p1_freq*p2_freq, (1-p1_freq)*(1-p2_freq))[negLD]
            corr_matrix/=tmp
            print("Saving LD",sname)
            with gzip.open(sample_location+'%s_LD_matrix.dat.gz'%sname, 'w') as fh:
                np.savetxt(fh, corr_matrix)

    for si,sname in enumerate(snames):
        with gzip.open(sample_location+'%s_LD_matrix.dat.gz'%sname, 'r') as fh:
            corr_matrix = np.loadtxt(fh)
        m = axs[si,1].matshow(corr_matrix, vmin=-1, vmax=1, cmap='RdBu_r')
        if si==2:
            plt.colorbar(m,cax=axs[-1,1], orientation='horizontal')
            axs[-1,1].set_xlabel('LD', fontsize=fs)

        axs[si,1].set_xlabel('SNV 1', fontsize=fs)
        axs[si,1].set_ylabel('SNV 2', fontsize=fs)
        axs[si,1].set_xticklabels([])
        axs[si,1].set_yticklabels([])


    plt.tight_layout()
    for fmt in ['.pdf', '.png']:
        plt.savefig('figs/dual_infections_TRIAL_Fig5'+fmt)



def plot_band_matrix(M, pos_map=None, upper=False):
    non_zero = np.where(M)
    coords = []
    for pi1, pi2 in zip(non_zero[0], non_zero[1]):
        p1 = pi1 if pos_map is None else pos_map[pi1]
        p2 = pi2 if pos_map is None else pos_map[pi2]
        if upper and p1>p2:
            continue
        px = 0.5*(p1+p2)
        py = (p2-p1)*np.sqrt(0.5)
        coords.append([px,py, M[pi1,pi2], p1,p2])

    coords = np.array(coords)

    from matplotlib.cm import RdBu_r as cmap
    #plt.scatter(coords[:,0], coords[:,1], c=coords[:,2])
    for px,py, v, p1, p2 in coords:
        plt.plot([p1,p2], [py,py], c=cmap((v+1)/2.0))

