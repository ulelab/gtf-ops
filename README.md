# gtf-ops

[![DOI](https://zenodo.org/badge/442399085.svg)](https://zenodo.org/badge/latestdoi/442399085)

**FiterGtf module**: Filtering GENCODE or ENSEMBLE annotation in GTF format for tag "basic" (first step) and "transcript_support_level" (second step).

🔴 **Important note: Currently, these tags are only available for Human (Hs) and Mouse (Mm) annotations.**

**ResolveUnnanotated module**: Annotates missing regions in iCount genomic segmentation as "genic_other".

## Features

gtf-ops package contains 2 functions, which can be run as scripts, that complement iCount segmentation.

**FilterGtf.py** filters GENCODE or ENSEMBL genomic annotation in GTF format. It can be used prior to running iCount segment to improve
genome-level segmentation by removing lower-confidence trancripts and favoring transcripts of full protein-coding genes over ncRNA, 
where they overlap.

**ResolveUnannotated.py** annotates genome segments that are not annotated by iCount segmentation. Missing annotations occur when a region 
overlaps with a gene in GTF annotation, but relevant transcripts were removed during filtering. Such region is not annotated as "intergenic", because it
overlaps a gene, nor can it be assigned any other region (5'UTR, 3'UTR, CDS, intron or ncRNA), due to lack of relevant transcripts.

Scripts can be run via command-line interface, more details are given in respective README files.
