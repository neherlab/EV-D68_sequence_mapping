# EV-D68_sequence_mapping

## Required Other Repositories
You will need to clone the [SVVC](https://github.com/neherlab/SVVC) repository _within_ this repository for all the code to work. You should then have a folder within this repository called `SVVC`.

## Intended data structure:
This code is intended to work where there is a 'sister' folder alongside called 'Enterovirus_D68_raw_data' which contains subfolders (example name 'MIC3108' is the default in this code) of sequencing runs, which in turn contain folders of paired-end `fastq.gz` reads, one folder per sample.

When called, specify the sequencing run subfolder within 'Enterovirus_D68_raw_data' and all these samples will be processed. The first step of the `trim` rule will access the reads in the sister folder, trim them, and output them in the folder `samples_by_labID`, with one folder per sample.
This is where all further processing of the sample will take place - the original sister folder is untouched.

### Sequence ID renaming
If desired, to maintain original sequence ID or lab ID secrecy when doing further processing (or to just have more sensible names), a symlink folder can be easily created so you can access samples by alternative names. Create a folder called `samples_by_strain` in the same directory as `samples_by_labID`. Create a tab-delimited two-column file with the lab ID in the first column and the new name in the second column. 

From within the `samples_by_strain` folder, use the command:
```
sed 's/^/ln -sf ..\/samples_by_labID\//g' < ../name_translation_table.tsv | sh
```

Such a tab-delimited file can also be used to create symlink folders with useful subsets of your data - for example, one sample per patient, or only sequences from a particular year (if only those are included in the tab-delimited file).

Here, the scripts that generate figures and analyses (in `analysis_scripts`) use the strain name of samples to access them via the symlink. This makes our code much easier to share, as it does not then contain any lab IDs, but instead the sample names we use in our paper and on GenBank.

### Mapping to Reference
We found we achieved much higher coverage across the genome when we mapped sequences to a reference within the same subclade. Thus, we do an initial mapping of 100 reads against references from the B1, B3, A1, and A2 subclades, and then map against the highest scoring reference. This is automatic within the pipeline. You can check the `references_performance.csv` file that's generated for each sample to see what reference was used. You can see the current references mapped against in the `references_M` folder.

