

import glob, re
# To run on cluster, call via nohup
# nohup snakemake --config 'run="801457"' --jobscript submit.sh --jobs 32 --cluster sbatch 1>log &
# If getting errors during mapping, try "export LMOD_DISABLE_SAME_NAME_AUTOSWAP="no"" to allow BWA and SAMtools to load at same time...

# To specify run folder, call as:
# snakemake --config 'run="801457"'

configfile: "snakemake_config.json"

regex = re.compile('Q-.|TA-.') # Do not process extra QC stuff
all_samples = [fname.split('/')[-1] for fname in glob.glob(config["raw_data"]+config["run"]+"/*") if not re.match(regex, fname.split('/')[-1])]
print(all_samples)

rule all:
    input:
        # trimmed = expand('samples_by_labID/{sample}/trimmed_reads_2.fq.gz', sample=all_samples)
        # ref = expand("samples_by_labID/{sample}/references_performance.csv", sample=all_samples)
        # mapped = expand('samples_by_labID/{sample}/mapped_reads.bam', sample=all_samples)
        # piled = expand("samples_by_labID/{sample}/" + "allele_counts.npz", sample=all_samples)
        consensed = expand("samples_by_labID/{sample}/consensus.fasta", sample=all_samples),
        figed = expand("samples_by_labID/{sample}/figures/coverage.png", sample=all_samples)
        

rule trim:
    output:
        "samples_by_labID/{sample}/trimmed_reads_1.fq.gz",
        "samples_by_labID/{sample}/trimmed_reads_2.fq.gz"
    params:
        min_length = 80,
        min_length_single = 90,
        run = config["run"],
        raw_data = config["raw_data"]
    shell:
        "ml Trim_Galore &&"
	    "trim_galore --length {params.min_length} --output samples_by_labID/{wildcards.sample} --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {params.raw_data}{params.run}/{wildcards.sample}/*_?.fastq.gz &&"
        "cat samples_by_labID/{wildcards.sample}/*val_1.fq.gz > {output[0]} &&"
        "cat samples_by_labID/{wildcards.sample}/*val_2.fq.gz > {output[1]} &&"
        "rm  samples_by_labID/{wildcards.sample}/*val_?.fq.gz"

#new rule for deciding on mapping reference to use
rule find_reference:
    input:
        rules.trim.output
    output:
        #references_performance.csv contains for each reference the error percentage achieved when mapping 1000 reads
        "samples_by_labID/{sample}/references_performance.csv"
    shell:
        "python ./references_M/reference_selection.py --trimmed_reads {input[0]} --ref_dir ./references_M/ --output {output}"

rule map:
    input:
        rules.find_reference.output,
        rules.trim.output
    output:
        "samples_by_labID/{sample}/mapped_reads.bam"
    shell:
        # Extracts the best reference from the find_reference output file
        """
        chosen_ref=$(python ./references_M/pick_top_reference.py {input[0]})
        echo $chosen_ref
        ml SAMtools &&
        ngm -r $chosen_ref -1 {input[1]} -2 {input[2]} -t 4 |samtools view -Sb - > {output}
        """
        # ml BWA SAMtools &
        #"bwa mem -k 6 {input} |samtools view -Sb - > {output}"
        #"minimap2 -ax sr {input} |samtools view -Sb - > {output}"

rule pileup:
    input:
        rules.map.output
    output:
        "samples_by_labID/{sample}/" + "allele_counts.npz" #%config["segment"]
    params:
        path_to_script = config["SVVC_dir"] + '/src',
        out_dir = "samples_by_labID/{sample}",
        primers = config["primers"] #"../../../primers.csv"
    shell:
        "python2 {params.path_to_script}/create_allele_counts.py --bam_file {input} --out_dir {params.out_dir} --primers {params.primers}"

#pair frequencies rule was here previously - but should only be run for duplicate samples

rule consensus:
    input:
        "samples_by_labID/{sample}/" + "allele_counts.npz" #%config["segment"]
    output:
        "samples_by_labID/{sample}/consensus.fasta",
        "samples_by_labID/{sample}/figures/coverage.png",
        "samples_by_labID/{sample}/figures/diversity.png",
        "samples_by_labID/{sample}/minor.fasta"
    params:
        path_to_script = config["SVVC_dir"] + '/src',
        out_dir = "samples_by_labID/{sample}",
        min_freq = 0.01,
        min_cov = 1000,
        all_counts = '--all_counts'
    shell:
        """
        echo {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} {params.all_counts} &
        python2 {params.path_to_script}/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} {params.all_counts} &
        echo {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir}  --min_freq {params.min_freq} --min_cov {params.min_cov} &
        python2 {params.path_to_script}/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} --min_freq {params.min_freq} --min_cov {params.min_cov}
        """
