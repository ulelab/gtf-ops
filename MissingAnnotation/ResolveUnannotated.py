import pandas as pd
import pybedtools as pbt
import csv
from plumbum.cmd import sort
import tempfile
import argparse

def cli():
    parser = argparse.ArgumentParser(description='Annotate genomic regions that are unnanotated in genomic iCount segmentation.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-s', '--segmentation', type=str, required=True,
                        help='iCount genome level segmentation in GTF format i.e "regions.gtf.gz".')
    required.add_argument('-a', '--annotation', type=str, required=True,
                        help='Annotation file from GENCODE or ENSEMBL in GTF format, which was as input for iCount segmentation.')
    required.add_argument('-fai', '--fasta_index', type=str, required=True,
                        help='Fasta index file generated with samtools.')
    required.add_argument('-o', '--outputdir', type=str, required=True,
                        help='Path to output folder.', default='.')
    args = parser.parse_args()
    print(args)
    return(args.segmentation, args.annotation, args.fasta_index, args.outputdir)


def ReadGtf(segmentation):
    df_segment = pd.read_csv(segmentation,
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
    return df_segment


def Fai2Bed(fai):
    df_chromosomes = pd.read_csv(fai, sep='\t', header=None, names=['chr', 'end', 'offset', 'linebases', 'linewidth'])
    df_chromosomes = df_chromosomes[['chr', 'end']].assign(start=0, name='.', score=0)
    df_chromosomes_p = df_chromosomes.copy()
    df_chromosomes_p['strand'] = '+'
    df_chromosomes_p = df_chromosomes_p[['chr', 'start', 'end', 'name', 'score', 'strand']]

    df_chromosomes_m = df_chromosomes.copy()
    df_chromosomes_m['strand'] = '-'
    df_chromosomes_m = df_chromosomes_m[['chr', 'start', 'end', 'name', 'score', 'strand']]

    df_chromosomes = pd.concat([df_chromosomes_p, df_chromosomes_m], ignore_index=True)
    bed_chr = pbt.BedTool.from_dataframe(df_chromosomes).sort()
    return(bed_chr)

def run(iCount_segment, gtf_annotation, fai, outputdir):
    # Read iCount genomic segment and convert it from GTF to BED format.
    print('Reading genomic segmentation.')
    df_segment = ReadGtf(iCount_segment)
    bed_segment = df_segment.assign(start=df_segment['start']-1, score=0)[['chrom', 'start', 'end', 'name', 'score','strand']]
    bed_segment = pbt.BedTool.from_dataframe(bed_segment).sort()
    # Convert fasta index to BED format - one entry spans one chromosome.
    bed_fai = Fai2Bed(fai)
    # Read annotation GTF, keep only gene-level entries and convert it to BED format.
    df_annotation = ReadGtf(gtf_annotation)
    df_annotation = df_annotation.loc[df_annotation['feature']=='gene']
    bed_annotation = df_annotation.assign(start=df_annotation['start']-1, score=0)[['chrom', 'start', 'end', 'annotations', 'score', 'strand']]
    bed_annotation = pbt.BedTool.from_dataframe(bed_annotation).sort()
    # Find regions that are unannotated in the iCount genome segmentation.
    print('Getting unannotated regions...')
    bed_missing = bed_fai.subtract(bed_segment, s=True).sort()
    print(f'Found {len(bed_missing)} unannotated genomic regions.')
    # Map gene annotation to the missing regions using bedtools and convert it to GTF format.
    print('Annotating regions...')
    bed_missing = bed_missing.map(bed_annotation, s=True, c=4, o='collapse')
    df_unnanotated = pd.read_csv(bed_missing.fn, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'annotations'])
    # Feature is genic_other.
    df_unnanotated = df_unnanotated.assign(feature='genic_other',
                                       start=df_unnanotated['start'] + 1,
                                       source='.',
                                       name2='.')
    df_unnanotated = df_unnanotated[['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations']]
    #Add missing regions to original iCount segment.
    print('Adding unnanotated regions to iCount segment as "genic_other".')
    df_segment = pd.concat([df_segment, df_unnanotated], ignore_index=True)
    # Sort GTF segment and write it to file
    outfile = f"{outputdir}/sorted.{iCount_segment.split('/')[-1].replace('.gz', '')}"
    with tempfile.NamedTemporaryFile(mode='w') as tmpfile:
        df_segment.to_csv(tmpfile.name, index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)
        cmd = (sort["-t\t", "-k1,1", "-k4,4n", tmpfile.name]) > outfile
        print(cmd())
    print(f'Saved the segment as {outfile}')
    return 0

def main():
    (iCount_segment, gtf_annotation, fai, outputdir) = cli()
    run(iCount_segment, gtf_annotation, fai, outputdir)

if __name__ == '__main__':
    main()
