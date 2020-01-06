#  Reference selection pipeline

This README files documents the python scripts used for incorporation of the automated reference selection into the mapping pipeline.

## reference_selection.py
This python script is called in order to evaluate the quality of mapping a small number of reads to all provided references.
### Input:
* **trimmed_reads**: Name of the file containing the trimmed reads used for mapping. Only one source is accepted, thus either the reads1 or reads2 have to be passed as argument, but not both.
* **ref_dir**: Directory where the references are stored
* **output**: Name of the file where the reference error rate should be stored
### Output:
* **trimmed_reads_1_top_1000.fq**: The first 1000 reads extracted from the trimmed_reads input file. Automatically placed into the same folder the **trimmed_reads** input is in.
* **mapped_reads_top_1000.sam**: The results after mapping the 1000 reads to the references. Will be evaluated directly after creation for each reference and is therefore overwritten by each new reference. Automatically placed into the same folder the **trimmed_reads** input is in.
* **references_performance**: CSV file listing for each reference file the error rates achieved during mapping of the 1000 reads, the references are sorted by performance. Stored as the file given in **output**.

## pick_top_reference.py
This python scripts takes the **references_performance** output from the script **reference_selection.py** and extracts the first reference (i.e. the best performing reference) from the sorted output file and returns this reference filename.
### Input:
* **references_performance**: CSV output from **reference_selection.py**
### Output:
* **reference_filename**: Returns the reference filename found with the lowest error rate.

