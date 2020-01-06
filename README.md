# EV-D68_sequence_mapping

## Intended data structure:
This code is intended to work where there is a 'sister' folder alongside called 'Enterovirus_D68_raw_data' which contains subfolders (example name 'MIC3108' is the default in this code) of sequencing runs, which in turn contain folders of paired-end `fastq.gz` reads, one folder per sample.

When called, specify the sequencing run subfolder within 'Enterovirus_D68_raw_data' and all these samples will be processed. The first step of the `trim` rule will access the reads in the sister folder, trim them, and output them in the folder 'samples_by_labID', with one folder per sample.
This is where all further processing of the sample will take place - the original sister folder is untouched.

If desired, to main original sequence or lab ID secrecy when doing further processing, a symlink folder can be easily created so you can access samples by alternative names.
