import glob, os, shutil
from Bio import SeqIO
import treetime

sample_location = 'samples_by_strain/'

#none of these contain QC samples
all_samples = glob.glob(sample_location+'EVD68*/consensus.fasta')
samples_14_16 = glob.glob(sample_location+'EVD68*_140*/consensus.fasta')
samples_14_16.extend( glob.glob(sample_location+'EVD68*_160*/consensus.fasta') )

#for including minor variants of dual infected:
duals = ["EVD68_SWE_045_160831_NFLG", "EVD68_SWE_046_160904_NFLG", "EVD68_SWE_007_140908_NFLG",
         "EVD68_BEL_009_181003_NFLG", "EVD68_CHE_006_181003_NFLG", "EVD68_ESP_004_180919_NFLG",
         "EVD68_NLD_001_180908_NFLG", "EVD68_BEL_004_180928_NFLG", "EVD68_NLD_007_181008_NFLG"]
all_samples_inc_dual = glob.glob(sample_location+'EVD68*/consensus.fasta')
for d in duals:
    all_samples_inc_dual.append(sample_location+d+'/minor.fasta')

all_samples_inc_dual_bkgrd = glob.glob(sample_location+'EVD68*/consensus.fasta')
for d in duals:
    #if d != "EVD68_SWE_007_140908_NFLG":
    all_samples_inc_dual_bkgrd.append(sample_location+d+'/minor.fasta')

all_samples_inc_all_minor_bkgrd = glob.glob(sample_location+'EVD68*/consensus.fasta')
all_samples_inc_all_minor_bkgrd.extend( glob.glob(sample_location+'EVD68*/minor.fasta') )

align_to_make = [('all_sequences', all_samples), ('sequences-14-16', samples_14_16),
                 ('all_seq_inc_dual', all_samples_inc_dual),
                 ('all_seq_inc_dual_bkgd', all_samples_inc_dual_bkgrd),
                 ('all_seq_all_minor_bkgd', all_samples_inc_all_minor_bkgrd)]

for outf, samps in align_to_make:
    outfname = 'alignments/'+outf+'.fasta'
    outfname_aln = 'alignments/'+outf+'_aln.fasta'
    outfname_tree = 'alignments/'+outf+'.nwk'

    with open(outfname, 'w') as ofile:
        for sample in samps:
            for seq in SeqIO.parse(sample, 'fasta'):
                ##REPLACE SEQ NAME
                seq.id = sample.split('/')[1]
                if 'minor' in sample:
                    seq.id = seq.id+'_minor'
                seq.description = ''
                SeqIO.write(seq, ofile, 'fasta')
        if outf == "all_seq_inc_dual_bkgd" or outf == "all_seq_all_minor_bkgd":
            for seq in SeqIO.parse("alignments/entero_background.fasta", 'fasta'):
                SeqIO.write(seq, ofile, 'fasta')

    os.system('mafft '+outfname + 	' > ' + outfname_aln)
    treetime.utils.tree_inference(outfname_aln, outfname_tree, methods=['fasttree'])

#make a copy with better name so next script is easier..
shutil.copyfile('alignments/all_seq_inc_dual_bkgd.nwk', 'alignments/all_seq_inc_dual_bkgd_aln_genome.nwk')

from Bio import AlignIO
aln = AlignIO.read('alignments/all_seq_inc_dual_bkgd_aln.fasta', 'fasta')


#make tree of F4 fragment
bound = (5302,7315)
ofile = "alignments/all_seq_inc_dual_bkgd_aln_f4"
AlignIO.write(aln[:,bound[0]:bound[1]], ofile+'.fasta', 'fasta')
treetime.utils.tree_inference(ofile+'.fasta', ofile+'.nwk', methods=['fasttree'])

#tree of F2 fragment
bound = (1812,3805)
ofile = "alignments/all_seq_inc_dual_bkgd_aln_f2"
AlignIO.write(aln[:,bound[0]:bound[1]], ofile+'.fasta', 'fasta')
treetime.utils.tree_inference(ofile+'.fasta', ofile+'.nwk', methods=['fasttree'])

#tree of F3 fragment
bound = (3571,5574)
ofile = "alignments/all_seq_inc_dual_bkgd_aln_f3"
AlignIO.write(aln[:,bound[0]:bound[1]], ofile+'.fasta', 'fasta')
treetime.utils.tree_inference(ofile+'.fasta', ofile+'.nwk', methods=['fasttree'])

#tree of F1 fragment
bound = (70,2210)
ofile = "alignments/all_seq_inc_dual_bkgd_aln_f1"
AlignIO.write(aln[:,bound[0]:bound[1]], ofile+'.fasta', 'fasta')
treetime.utils.tree_inference(ofile+'.fasta', ofile+'.nwk', methods=['fasttree'])

os.remove("tmp.nwk")
os.remove("fasttree_stderr")
