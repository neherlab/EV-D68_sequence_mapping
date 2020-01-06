import subprocess
import pysam
import numpy as np
from Bio import SeqIO
import glob, sys,os, argparse

if __name__ == '__main__':
    import pandas as pd

    # Parse input args
    parser = argparse.ArgumentParser(description='map a few reads to every reference and evaluate quality',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--trimmed_reads', required=True, type=str, help='the file containing the trimmed reads')
    parser.add_argument('--ref_dir', required=True, type=str, help='directory to references')
    parser.add_argument('--output', required=True, type=str, help='filename of reference quality output')

    args = parser.parse_args()

    n = 1000 # Number of reads to map
    ref_dir = args.ref_dir
    references = [ref_dir + file for file in os.listdir(ref_dir) if file.endswith(".fasta")] # List of reference files
    # .fasta ending required since ngm indexes the references

    trimmed_reads = args.trimmed_reads
    data_path = trimmed_reads[:(trimmed_reads.rfind('/') + 1)] # Folder containig the trimmed read will be used for various output

    percentage_of_mismatches = {} # Dictionary containing for every reference the error rate
    for reference in references:
        with open(data_path + "trimmed_reads_1_top_" + str(n) + ".fq", "w") as output_file: # Extract top n reads from the trimmed reads file
            trimmed_reads1_all = subprocess.Popen(["gzip", "-cd", trimmed_reads], stdout=subprocess.PIPE)
            trimmed_reads1 = subprocess.Popen(["head", "-" + str(n)], stdin=trimmed_reads1_all.stdout, stdout=output_file)
            trimmed_reads1_all.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        trimmed_reads1.wait() # Wait for this process to finish before continuing

        with open(data_path + "mapped_reads_top_" + str(n) + ".sam", "w") as output_file: # Map the n reads to the reference using ngm and store the result in the same folder
            mapping = subprocess.Popen(["ngm", "-q", data_path + "trimmed_reads_1_top_" + str(n) + ".fq", "-r", reference], stdin=trimmed_reads1.stdout, stdout=output_file)
        mapping.wait() # Wait for this process to finish before continuing

        # Compare mapped reads to the reference and extract error rates
        reference_seq = SeqIO.read(reference, "fasta").seq
        with pysam.Samfile(data_path + "mapped_reads_top_" + str(n) + ".sam") as samfile:
            percentage = []
            for read in samfile:
                query = read.seq
                mismatches = 0
                for (q,r) in read.get_aligned_pairs(matches_only=True):
                    if query[q] != reference_seq[r]:
                        mismatches += 1
                percentage.append(mismatches/len(query))

            percentage_of_mismatches[reference] = np.mean(percentage)

    # Store error rates in csv output file
    with open(args.output, "w") as out:
        for key in sorted(percentage_of_mismatches, key=percentage_of_mismatches.get):
            out.write(key + "\t" + str(percentage_of_mismatches[key]) + "\n")
