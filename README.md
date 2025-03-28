# Makassar_project

The scripts are a part of the project that analyzes the off-targets in Makassar locus.
1. The python scripts were used to process the fastq files, align them to the off-target sequences and calculate the frequency of all the bases at each position in the sequence.
   The base counts for each position will be used to calculate the per base error rate.

2. The R scripts were used to process the alignment results to calculate the per base error rate in the protospacer sequence,
   the background error rate for each animal based on the pre-treatment samples and to indentify substitutions of interest in the edited samples.
   Sequence to run through R scripts:
   1. Base error rate.R
   2. Background error rate.R
   3. Identify mutations.R
