import pandas as pd
import csv
import argparse

def cli():
    parser = argparse.ArgumentParser(description='Filter genomic annotation from GENCODE or ENSEMBL in GTF format.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-a', '--annotation', type=str, required=True,
                        help='Annotation file from GENCODE or ENSEMBL in GTF format.')
    required.add_argument('-o', '--outputdir', type=str, required=True,
                        help='Path to output folder.', default='.')
    args = parser.parse_args()
    print(args)
    return(args.annotation, args.outputdir)

def filter_gff(gtf_file, outputdir):
    # Parse gtf file into pandas dataframe.
    print("Reading annotation file.")
    annotation = pd.read_csv(gtf_file,
                             sep='\t',
                             names=['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations'],
                             header=None,
                             comment='#',
                             dtype={
                                    "chrom": str,
                                    "source": str,
                                    "feature": str,
                                    "start": int,
                                    "end": int,
                                    "name": str,
                                    "strand": str,
                                    "name2": str,
                                    "annotations": str,
                                    })
    # Filter transcripts by basic tag
    print("Number of entries in input annotation:", len(annotation))
    annotation = annotation.loc[annotation['annotations'].str.contains('tag "basic"') | (annotation['feature'] == 'gene'), :]
    print("Number of entries after filtering for tag \"basic\":", len(annotation))
    # Filter annotation gene-by-gene by transcript level support (TSL), to keep higher confidence transcripts where possible (TSL1 and 2)
    df_TSL = annotation.loc[annotation['annotations'].str.contains('transcript_support_level "1|transcript_support_level "2', regex=True), :]
    gene_ids = df_TSL["annotations"].str.split(";", n=1, expand=True)[0].unique().tolist()
    print("Number of genes that contain TSL1 or TSL2 transcripts:", len(gene_ids))
    # Keeping only TSL1 and TSL2 entries for genes that contain them, discardig other entries (no TSL information or TSL3-5)
    print('Filtering out low-confidence transcripts.')
    df_t = annotation.loc[(annotation['feature'] != 'gene') & (annotation['annotations'].str.contains('|'.join(gene_ids))) , :]
    df_t = df_t.loc[~df_t['annotations'].str.contains('transcript_support_level "1"|transcript_support_level "2"', regex=True)]
    annotation.drop(index=df_t.index, inplace=True)
    print("Number of entries in filtered annotation.", len(annotation))
    print('Saving filtered gtf file.')
    annotation.to_csv(f"{outputdir}/filtered.{gtf_file.split('/')[-1]}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)

def main():
    (gtf_file, outputdir) = cli()
    filter_gff(gtf_file, outputdir)

if __name__ == '__main__':
    main()