import pandas as pd
import csv
import argparse

def cli():
    parser = argparse.ArgumentParser(description='Filter genomic annotation from GENCODE or ENSEMBL in GTF format by tag \"basic\" and transcript_support_level.'\
    ' These flags are currently (20250312) annotated for Homo sapiens and Mus musculus organisms.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-a', '--annotation', type=str, required=True,
                        help='Annotation file from GENCODE or ENSEMBL in GTF format.')
    required.add_argument('-o', '--outputdir', type=str, required=True,
                        help='Path to output folder.')
    args = parser.parse_args()
    print(args)
    return(args.annotation, args.outputdir)

def filter_gff(gtf_file, outputdir):
    # Parse gtf file into pandas dataframe.
    print("Reading annotation file.")
    input_annotation = pd.read_csv(gtf_file,
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
    # Get gene ids from input annotation
    input_gene_ids = set(input_annotation.loc[input_annotation['feature'] == 'gene', 'annotations'].str.split(";", n=1, expand=True)[0].unique().tolist())
    # Filter transcripts by basic tag
    print("Number of entries in input annotation:", len(input_annotation))
    # Check if annotation contains tag "basic"
    print("Checking for basic flag...")
    # Made basic more specific in case it would be part of gene-name or other annotation
    basic = input_annotation['annotations'].str.contains('tag "basic"', regex=True)
    if basic.any():
        print("Basic flag available.")
        nbasic = basic.value_counts()[True]
        print(f"{nbasic} entries flagged as basic.")
        # Filter annotation for tag "basic", but keep all gene entries
        annotation = input_annotation.loc[input_annotation['annotations'].str.contains('tag "basic"') | (input_annotation['feature'] == 'gene'), :].copy()

        # Check if every gene contains at least one transcript entry
        print("Checking if every gene contains at least one transcript entry: each gene should have a basic set of transcripts.")
        gene_ids_transcripts = set(annotation.loc[annotation['feature'] == 'transcript', 'annotations'].str.split(";", n=1, expand=True)[0].unique().tolist())
        # Print the first five gene ids for visual
        print("First five gene ids:", list(gene_ids_transcripts)[:5])
        if input_gene_ids == gene_ids_transcripts:
            print("Every gene contains at least one transcript entry.")
        else:
            print("Some genes do not contain any transcript entry tagged as basic - this indicates that the input gtf may not be correctly annotated or that some entries are missing.")
            print('Returning input annotation as output. Exiting.')
            input_annotation.to_csv(f"{outputdir}/unfiltered.{gtf_file.split('/')[-1]}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
            return

        print("Number of entries after filtering for tag \"basic\":", len(annotation))

        # Filter annotation gene-by-gene by transcript level support (TSL), to keep higher confidence transcripts where possible (TSL1 and 2)
        print("Checking for transcript support level flag...")
        if annotation['annotations'].str.contains('transcript_support_level', regex=True).any():
            print("transcript_support_level flag available.")
            if annotation['annotations'].str.contains('transcript_support_level "1|transcript_support_level "2', regex=True).any():
                print("TSL1 or TSL2 flag available.")
                # Get only TSL1 and TSL2 entries to perform elimination on
                df_TSL = annotation.loc[annotation['annotations'].str.contains('transcript_support_level "1|transcript_support_level "2', regex=True), :]
                # Get gene ids that contain TSL1 or TSL2 transcripts
                gene_ids_highconf = df_TSL["annotations"].str.split(";", n=1, expand=True)[0].unique().tolist()
                print("First five gene ids that contain TSL1 or TSL2 transcripts:", gene_ids_highconf[:5])
                
                print("Number of genes that contain TSL1 or TSL2 transcripts:", len(gene_ids_highconf))
                # Keeping only TSL1 and TSL2 entries for genes that contain them, discardig other entries (no TSL information or TSL3-5)
                print('Filtering out low-confidence transcripts.')
                # Get all entries that are not gene entries and contain gene ids that contain TSL1 or TSL2 transcripts
                df_t = annotation.loc[(annotation['feature'] != 'gene') & (annotation['annotations'].str.contains('|'.join(gene_ids_highconf))) , :]
                # Remove TSL1 and TSL2 entries, this leaves lower-confidence entries
                df_t = df_t.loc[~df_t['annotations'].str.contains('transcript_support_level "1|transcript_support_level "2', regex=True)]
                # Drop non-TSL1 / 2 from annotation
                annotation.drop(index=df_t.index, inplace=True)
                # Check that each gene contains at least one transcript entry
                print("Checking if every gene contains at least one transcript entry after TSL filtering.")
                # Get gene ids of retained transcripts
                filtered_gene_ids = set(annotation.loc[annotation['feature'] == 'transcript', 'annotations'].str.split(";", n=1, expand=True)[0].unique().tolist())
                if input_gene_ids == filtered_gene_ids:
                    print("Every gene contains at least one transcript entry.")
                else:
                    print("Some genes do not contain any transcript entry - this indicates that the input gtf may not be correctly annotated or that some entries are missing.")
                    print('Returning input annotation as output. Exiting.')
                    input_annotation.to_csv(f"{outputdir}/unfiltered.{gtf_file.split('/')[-1]}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
                    return

            else:
                print("No TSL1 or TSL2 flag available.")
                print("Keeping all transcripts flagged as basic.")
        else:
            print("No transcript_support_level flag available.")
            print("Keeping all transcripts flagged as basic.")
        
        print("Number of entries in filtered annotation.", len(annotation))
        print('Saving filtered gtf file.')
        annotation.to_csv(f"{outputdir}/filtered.{gtf_file.split('/')[-1]}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
        return
    else:
        print('No tag \"basic\". Returning input annotation as output. Exiting.')
        input_annotation.to_csv(f"{outputdir}/unfiltered.{gtf_file.split('/')[-1]}", header=None, index=None, sep='\t', quoting=csv.QUOTE_NONE)
        return

def main():
    (gtf_file, outputdir) = cli()
    filter_gff(gtf_file, outputdir)

if __name__ == '__main__':
    main()
