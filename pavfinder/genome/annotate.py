from pybedtools import BedTool, create_interval_from_list
from itertools import groupby
from intspan import intspan
import multiprocessing as mp
import sys


def get_acen_coords(cytobands_file):
    """Extracts acentromere coordinates from UCSC cytobands file

    column 5 = "acen"

    Args:
        cytobands_file: (str) path to cytobands file
    Returns:
        Dictionary where key = chrom
                         value = list of (start, end)
    """
    acen_coords = None
    try:
        fh = open(cytobands_file, 'r')
        with fh:
            acen_coords = {}
            for line in fh:
                chrom, start, end, band, stain = line.rstrip('\n').split('\t')
                if stain == 'acen':
                    try:
                        acen_coords[chrom].append((start, end))
                    except BaseException:
                        acen_coords[chrom] = [(start, end)]

    except BaseException:
        sys.stderr.write("can't open cytobands file:{}\n".format(cytobands_file))

    return acen_coords


def update_features(variant, features):
    """Updating variant's annotation given features

    Args:
        variant: (Variant) Variant to update annotation
        features: (tuple) 2 Pybetools Interval objects that overlap the 2 breakpoints
    """
    for i in range(2):
        if features[i] is not None:
            variant.genes[i] = features[i].attrs['gene_name']
            variant.transcripts[i] = features[i].attrs['transcript_id']
            variant.gene_strands[i] = features[i].strand
            if 'exon_number' in features[i].attrs:
                variant.exons[i] = int(features[i].attrs['exon_number'])
            if 'intron_number' in features[i].attrs:
                variant.introns[i] = int(features[i].attrs['intron_number'])
            if 'exon_bound' in features[i].attrs:
                if features[i].attrs['exon_bound'] == 'True':
                    variant.exon_bounds[i] = True
                else:
                    variant.exon_bounds[i] = False


def annotate_rna_event(variant):
    if variant.genes[0] != 'NA' and variant.genes[1] != 'NA':
        if variant.genes[0] != variant.genes[1]:
            variant.rna_event = 'gene_fusion'
        elif variant.exons[0] != 'NA' and variant.exons[0] == variant.exons[1] and variant.rearrangement == 'dup':
            variant.rna_event = 'ITD'
        elif variant.exons[0] != 'NA' and variant.exons[0] != variant.exons[1] and\
                variant.exon_bounds[0] and variant.exon_bounds[1] and\
                (variant.rearrangement == 'dup' or variant.rearrangement == 'trp'):
            variant.rna_event = 'PTD'


def annotate_gene_fusion(variant):
    if variant.orients[0] == 'L' and variant.orients[1] == 'L':
        if variant.gene_strands[0] == '+' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'antisense'
        elif variant.gene_strands[0] == '+' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[0], variant.genes[1]
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[1], variant.genes[0]
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'antisense'

    elif variant.orients[0] == 'R' and variant.orients[1] == 'R':
        if variant.gene_strands[0] == '+' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'antisense'
        elif variant.gene_strands[0] == '+' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[1], variant.genes[0]
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[0], variant.genes[1]
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'antisense'

    elif variant.orients[0] == 'L' and variant.orients[1] == 'R':
        if variant.gene_strands[0] == '+' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[0], variant.genes[1]
        elif variant.gene_strands[0] == '+' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'antisense'
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'antisense'
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[1], variant.genes[0]

    elif variant.orients[0] == 'R' and variant.orients[1] == 'L':
        if variant.gene_strands[0] == '+' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[1], variant.genes[0]
        elif variant.gene_strands[0] == '+' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'antisense'
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '+':
            variant.fusion_type = 'antisense'
        elif variant.gene_strands[0] == '-' and variant.gene_strands[1] == '-':
            variant.fusion_type = 'sense'
            variant.gene5, variant.gene3 = variant.genes[0], variant.genes[1]

    if variant.fusion_type == 'sense':
        if variant.exon_bounds[0] == variant.exon_bounds[1] and variant.exon_bounds[0]:
            variant.frames = 'in'


def overlap(bed_file, gtf_file, result_file):
    variants_bed = BedTool(bed_file)
    gtf = BedTool(gtf_file)

    variants_bed.intersect(gtf, wb=True).moveto(result_file)


def pick_feature(break_pt, orient, features):
    """Pick the 'best' feature for a given breakpoint given a list of features
    This is the workflow to determine the 'best':
    is it(breakpoint) in an exon?
    if yes, does it coincide with an exon boundary
    if yes, is the transcript coding
    if it's not in an exon, is it in an intron
    if yes, is the trancript coding?
    if it's not an exon or intron, returns None
    if more than 1 possible choices,
    sort by transcript_id and then exon_number, returns the first

    Args:
        break_pt: (tuple) (chromosome, position) of breakpoint
        orient: (str) 'L' or 'R'
        features: (list) Features that overlap breakpoint

    Returns:
        'best' feature (Interval object of pyBedTools)
    """
    best_features = []
    exons = [feature for feature in features if feature[2] == 'exon']
    if exons:
        for feature in exons:
            feature.attrs['exon_bound'] = 'False'

        exon_bounds = [exon for exon in exons if (break_pt[1] == exon.start + 1 and orient == 'R') or
                                                 (break_pt[1] == exon.end and orient == 'L')]

        if exon_bounds:
            for feature in exon_bounds:
                feature.attrs['exon_bound'] = 'True'

            coding_exons = [exon for exon in exon_bounds if 'gene_biotype' in exon.attrs and
                            exon.attrs['gene_biotype'] == 'protein_coding']

            if coding_exons:
                best_features = coding_exons
            else:
                best_features = exon_bounds

        else:
            best_features = exons

        if len(best_features) > 1:
            by_transcript_id = sorted(best_features, key=lambda feature: feature['transcript_id'])
            by_exon_num = sorted(by_transcript_id, key=lambda feature: int(feature['exon_number']), reverse=True)
            return by_exon_num[0]

        elif len(best_features) == 1:
            return best_features[0]

    else:
        introns = [feature for feature in features if feature[2] == 'intron']

        if introns:
            coding_introns = [intron for intron in introns if 'gene_biotype' in intron.attrs and
                              intron.attrs['gene_biotype'] == 'protein_coding']

            if coding_introns:
                best_features = coding_introns
            else:
                best_features = introns

        if len(best_features) > 1:
            by_transcript_id = sorted(best_features, key=lambda feature: feature['transcript_id'])
            by_intron_num = sorted(by_transcript_id, key=lambda feature: int(feature['intron_number']), reverse=True)
            return by_intron_num[0]

        elif len(best_features) == 1:
            return best_features[0]

    return None


def locate_features(breaks, orients, features):
    """Find the 'best' gene features of the breakpoints of a given event
    It will first determine which features overlap which breakpoint.
    If there are features that overlap both breakpoints, only they will be considered
    in picking the 'best' suitable one.

    Args:
        breaks: (tuple) the 2 breakpoints ((chr1, pos1), (chr2, pos2))
        orients: (tuple) the 2 orientations ('L|R', 'L|R')
        features: (list) Interval objects of a given event from parsing the bedpe overlap file

    Returns:
        A tuple of 2 feature (interval object) picked to annotate the 2 breakpoints
        can be (None, None) if nothing is found)
    """
    # use intspan to intersect individual breakpoint with feature coodinates
    break1_span = intspan('{}-{}'.format(breaks[0][1], breaks[0][1]))
    break2_span = intspan('{}-{}'.format(breaks[1][1], breaks[1][1]))

    # categories features where 1, 2, or both breakpoints overlap
    # use Set because there may be redundancy
    overlaps = {'both': set(), '1': set(), '2': set()}

    for feature in features:
        feature_span = intspan('{}-{}'.format(feature.start + 1, feature.stop))

        overlap1 = True if feature.chrom == breaks[0][0] and feature_span & break1_span else False
        overlap2 = True if feature.chrom == breaks[1][0] and feature_span & break2_span else False
        if overlap1 and overlap2:
            overlaps['both'].add(feature)
        elif overlap1:
            overlaps['1'].add(feature)
        elif overlap2:
            overlaps['2'].add(feature)

    # only considers features that overlap both breakpoints if such are found
    if overlaps['both']:
        best_feature1 = pick_feature(breaks[0], orients[0], overlaps['both'])
        best_feature2 = pick_feature(breaks[1], orients[1], overlaps['both'])
    else:
        best_feature1 = pick_feature(breaks[0], orients[0], overlaps['1'])
        best_feature2 = pick_feature(breaks[1], orients[1], overlaps['2'])

    return best_feature1, best_feature2


def parse_overlaps(bed_file):
    """Parses bedtools results of event breakpoints and annotation file and returns single
    best feature for each breakpoint

    Args:
        bedpe_file: (str) Path of bedtools result of bedpe file of event coords vs. annotation

    Returns:
        A dictionary where key = name of variant as appeared in bedpe file
                           value = (feature1, feature2),
                                   where feature1, feature2 = features overlapping the 2 breakpoints
                                         and they are be None
    """
    count = 1
    results = {}
    with open(bed_file, 'r') as f:
        for variant_name, group in groupby(
                f, lambda line: line.split('\t')[6]):
            count += 1
            group = list(group)
            cols = group[0].split('\t')
            # break = (chrom, pos)
            break1 = cols[0], int(cols[2])
            break2 = cols[3], int(cols[5])

            orient1, orient2 = variant_name.split('-')[-2:]

            # extract all the overlapping features (Interval objects)
            features = []
            for line in group:
                cols = line.split('\t')
                features.append(create_interval_from_list(cols[-9:]))

            # pick the best 2 features
            if features:
                results[variant_name] = locate_features(
                    (break1, break2), (orient1, orient2), features)
            else:
                results[variant_name] = None, None
    f.close()

    return results


def parse_overlaps(bed_file, variant_keys=None):
    """Parses bedtools results of event breakpoints and annotation file and returns best features for
    each breakpoint

    Args:
        bed_file: (str) Path of bedtools results of event coords vs. annotation
        variant_keys: (Set) Keys of variants (variant.key()) to consider only

    Returns:
        List of (variant_key, (feature1, feature2))
        where a 'feature' is an Interval object (from Pybedtools) that's chosen for the breakpoint in question
    """
    results = []

    with open(bed_file, 'r') as f:
        for variant_key, group in groupby(f, lambda line: line.split('\t')[6]):
            if variant_keys is not None and variant_key not in variant_keys:
                continue

            group = list(group)
            cols = group[0].split('\t')
            # break = (chrom, pos)
            break1 = cols[0], int(cols[2])
            break2 = cols[3], int(cols[5])

            orient1, orient2 = variant_key.split('-')[-2:]

            # extract all the overlapping features (Interval objects)
            features = []
            for line in group:
                cols = line.split('\t')
                features.append(create_interval_from_list(cols[-9:]))

            result = None, None
            if features:
                result = locate_features((break1, break2), (orient1, orient2), features)

            results.append((variant_key, result))
    f.close()

    return results


def create_batches(bed_file, variants, size):
    """Creates batches for parallel parsing of bedtools results

    Args:
        bed_file: (str) Path of bedtools results of event coords vs. annotation
        variants: (list) Entire list of variant keys to be divided
        size: (int) Number of variant keys to be processed in each batch

    Yields:
        Tuple of (bed_file, list of variant_keys)
    """
    if size == 0:
        yield bed_file, variants
    else:
        for i in range(0, len(variants), size):
            if len(variants) - (i + size) < size:
                yield bed_file, variants[i:]
                break
            else:
                yield bed_file, variants[i:i + size]


def worker(args):
    """Creates worker process of multi-processing

    Args:
        args: (tuple) Bed file, list of variant keys

    Returns:
        Function call to parse_overlaps()
    """
    bed_file, variant_keys = args
    return parse_overlaps(bed_file, set(variant_keys))


def parallel_parse_overlaps(bed_file, variant_keys, num_procs):
    """Parses annotation overlap results using Multi-processing module

    Args:
        bed_file: (str) Path of bedtools results of event coords vs. annotation
        variant_keys: (list) Entire list of variant keys to be divided
        num_procs: (int) Number of processes the parallelization is going to use

    Returns:
        A dictionary of results with the variant key as the key
        and a tuple of Interval objects as the value
    """
    batches = list(create_batches(bed_file, variant_keys, len(variant_keys) / num_procs))
    pool = mp.Pool(processes=num_procs)
    batch_results = pool.map(worker, batches)
    pool.close()
    pool.join()

    results = {}
    for batch_result in batch_results:
        for variant_key, result in batch_result:
            results[variant_key] = result

    return results


def overlap_pe(variants_bedpe_file, gtf_file, result_file):
    variants_bed = BedTool(variants_bedpe_file)
    gtf = BedTool(gtf_file)
    olaps = variants_bed.pair_to_bed(gtf, stream=True).moveto(result_file)
