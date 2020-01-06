import glob, os
from Bio import SeqIO

sample_location = 'samples_by_strain2/'

#none of these contain QC samples
all_samples = glob.glob(sample_location+'EVD68*/consensus.fasta')
samples_18 = glob.glob(sample_location+'EVD68*_18*/consensus.fasta')

align_to_make = [('all_sequences_2018', all_samples), ('sequences-18', samples_18)]

for outf, samps in align_to_make:
    outfname = 'alignments/'+outf+'.fasta'
    outfname_aln = 'alignments/'+outf+'_aln.fasta'

    with open(outfname, 'w') as ofile:
        for sample in samps:
            for seq in SeqIO.parse(sample, 'fasta'):
                ##REPLACE SEQ NAME
                seq.id = sample.split('/')[1]
                seq.description = ''
                SeqIO.write(seq, ofile, 'fasta')

    os.system('mafft '+outfname + 	' > ' + outfname_aln)
