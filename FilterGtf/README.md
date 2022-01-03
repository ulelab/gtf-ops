# FilterGtf

FilterGtf.py filters GENCODE or ENSEMBL genomic annotation in GTF format. 

**Motivation**

Pre-filtering genomic annotation is crucial to ensure that accurate genome-level segmentation is produced by iCount segment. When iCount segmentation is applied to unfiltered annotation, this can result in introns of protein-coding genes to be annotated as ncRNA and in shifting of CDSs and UTRs to include low-confidence transcripts.

**Filtering rules**

1. No filtering is done on gene level.
2. Entries below gene-level tagged as "basic" are retained. The transcripts tagged as "basic" form part of a subset of representative transcripts for each gene. This subset prioritizes full-length protein-coding transcripts over partial or non-protein-coding transcripts within the same gene, and intends to highlight those transcripts that will be useful to the majority of users.
3. If a gene contains entries with transcript_support_level 1 or 2, only those transcripts are retained by filtering and other less-supported trascripts are discarded.

**Dependencies**:
```
python=3.7
pandas

Vesions used in development and testing:
pandas=0.24.2
```

**Set up**:

We recommend running FilterGtf in a conda environment so all the dependencies are managed for you, to set this up run the following command from your gtf-ops directory:
```
conda env create -f environment.yml
```
Before you run FilterGtf, activate your gtf-ops nvironment:
```
conda activate gtf-ops
```

**Usage**
```
python3 FilterGtf.py [-h] -a ANNOTATION -o OUTPUTDIR
Required arguments:
  -a, --annotation ANNOTATION: Annotation file from GENCODE or ENSEMBL in GTF format.
  -o, --outputdir OUTPUTDIR: Path to output folder.
```

**Outputs**

Filtered gtf file.

