import pandas as pd
import pybedtools as pbt
import csv
from plumbum.cmd import sort
import tempfile
import argparse

def cli():
    parser = argparse.ArgumentParser(description='Annotate genomic regions that are unnanotated in genomic iCount segmentation.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-fseg', '--filtered_segmentation', type=str, required=True,
                        help='iCount genome level segmentation in GTF format i.e "regions.gtf.gz", produced from filtered genomic annotation.')
    required.add_argument('-useg', '--unfiltered_segmentation', type=str, required=True,
                        help='iCount genome level segmentation in GTF format i.e "regions.gtf.gz", produced from raw genomic annotation.')
    required.add_argument('-a', '--annotation', type=str, required=True,
                        help='Raw annotation file from GENCODE or ENSEMBL in GTF format.')
    required.add_argument('-fai', '--fasta_index', type=str, required=True,
                        help='Fasta index file generated with samtools.')
    required.add_argument('-o', '--outputdir', type=str, required=True,
                        help='Path to output folder.', default='.')
    parser.add_argument('-go', '--genic_other', action='store_true')
    args = parser.parse_args()
    print(args)
    return(args.filtered_segmentation, args.unfiltered_segmentation, args.annotation, args.fasta_index, args.outputdir, args.genic_other)


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

def run(filtered_segment, unfiltered_segment, gtf_annotation, fai, outputdir, genic_other=False):
    # Read filtered iCount genomic segment and convert it from GTF to BED format.
    print('Reading genomic segmentation.')
    df_segment = ReadGtf(filtered_segment)
    bed_segment = df_segment.assign(start=df_segment['start']-1, score=0)[['chrom', 'start', 'end', 'feature', 'score','strand']]
    bed_segment = pbt.BedTool.from_dataframe(bed_segment).sort()
    # Read unfiltered iCount genomic segment and convert it from GTF to BED format.
    df_unfiltered = ReadGtf(unfiltered_segment)
    bed_unfiltered = df_unfiltered.assign(start=df_unfiltered['start']-1, score=0)[['chrom', 'start', 'end', 'feature', 'score','strand', 'annotations']]
    bed_unfiltered = pbt.BedTool.from_dataframe(bed_unfiltered).sort()

    # Convert fasta index to BED format - one entry spans one chromosome.
    bed_fai = Fai2Bed(fai)

    # Read annotation GTF, keep only gene-level entries and convert it to BED format.
    print('Getting gene-level annotation...')
    df_annotation = ReadGtf(gtf_annotation)
    df_annotation = df_annotation.loc[df_annotation['feature']=='gene']
    bed_annotation = df_annotation.assign(start=df_annotation['start']-1, score=0)[['chrom', 'start', 'end', 'annotations', 'score', 'strand']]
    bed_annotation = pbt.BedTool.from_dataframe(bed_annotation).sort()

    # Find regions that are unannotated in the iCount genome segmentation.
    print('Getting unannotated regions...')
    bed_missing = bed_fai.subtract(bed_segment, s=True,  nonamecheck=True).sort()
    print(f'Found {len(bed_missing)} unannotated genomic regions.')
    # Use intersect to split unnanotated regions
    intersect = bed_missing.intersect(bed_unfiltered, s=True, nonamecheck=True).sort()
    print('Annotating regions with gene information...')
    if not genic_other:
        # Intersect missing regions with unfiltered segment to get transcript region
        print('Annotating missing regions in iCount segment with transcript regions...')
        # Annotate with annotations (column 7) and feature (column 4)
        missingAnnotated = intersect.map(bed_unfiltered, s=True, c=[7,4], o='collapse', nonamecheck=True).sort()
        df_unnanotated = pd.read_csv(missingAnnotated.fn, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'annotations', 'feature'])
        df_unnanotated = df_unnanotated.assign(start=df_unnanotated['start'] + 1,
                                               source='.',
                                               name2='.')
    else:
        print('Annotationg missing regions in iCount segment as "genic_other".')
        # Annotate with annotations (column 7), feature is genic_other
        missingAnnotated = intersect.map(bed_unfiltered, s=True, c=7, o='collapse', nonamecheck=True).sort()
        df_unnanotated = pd.read_csv(missingAnnotated.fn, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'annotations'])
        # Feature is genic_other.
        df_unnanotated = df_unnanotated.assign(feature='genic_other',
                                           start=df_unnanotated['start'] + 1,
                                           source='.',
                                           name2='.')
        # Find regions that are missing from main annotation and annotate them with gene entries from annotations gtf file
        # Get complement of raw segment to find missing regions
        regsMissingFromMain = bed_fai.subtract(bed_unfiltered, s=True,  nonamecheck=True).sort()
        # Annotate them from annotation gtf (gene level annotation) and format entries to replace "gene_type" with "biotype"
        regsMissingFromMain = regsMissingFromMain.map(bed_annotation, s=True, c=4, o='collapse', nonamecheck=True)
        dfMissingFromMain = pd.read_csv(regsMissingFromMain.fn, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'annotations'])
        dfMissingFromMain['annotations'] = dfMissingFromMain['annotations'].apply(lambda x: str(x).replace('gene_type', 'biotype'))
        dfMissingFromMain = dfMissingFromMain.assign(feature='genic_other',
                                           start=dfMissingFromMain['start'] + 1,
                                           source='.',
                                           name2='.')
        dfMissingFromMain = dfMissingFromMain[['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations']]
        # Combine all missing regions
        df_unnanotated = pd.concat([df_unnanotated, dfMissingFromMain])
    df_unnanotated = df_unnanotated[['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations']]
    #Add missing regions to original iCount segment.
    print('Adding annotated missisng regions to iCount segment...')
    df_segment = pd.concat([df_segment, df_unnanotated], ignore_index=True)
    print('N segment entries:', len(df_segment))
    # Sort GTF segment and write it to file
    if genic_other:
        identifier = 'genic_other'
    else:
        identifier = 'annotated'
    outfile = f"{outputdir}/sorted.{identifier}.{filtered_segment.split('/')[-1].replace('.gz', '')}"
    with tempfile.NamedTemporaryFile(mode='w') as tmpfile:
        df_segment.to_csv(tmpfile.name, index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)
        cmd = (sort["-t\t", "-k1,1", "-k4,4n", tmpfile.name]) > outfile
        print(cmd())
    print(f'Saved the segment as {outfile}')
    return 0

def main():
    (filtered_segmentation, unfiltered_segmentation, gtf_annotation, fai, outputdir, genic_other) = cli()
    run(filtered_segmentation, unfiltered_segmentation, gtf_annotation, fai, outputdir, genic_other)

if __name__ == '__main__':
    main()
