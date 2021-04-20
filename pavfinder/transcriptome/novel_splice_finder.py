from .alignment import reverse_complement
from .adjacency import Adjacency
from .transcript import Transcript
from collections import OrderedDict
from copy import deepcopy

report_items = OrderedDict(
    [('event', 'event'),
     ('gene', None),
     ('transcript', None),
     ('size', 'size'),
     ('exon1', 'exons,0',),
     ('exon2', 'exons,1',),
     ('seq_id', 'seq_id'),
     ('seq_breaks', 'seq_breaks'),
     ('in_frame', 'in_frame'),
     ('splice_motif', 'splice_motif'),
     ('probe', 'probe'),
     ('support_reads', 'spanning'),
     ('splice_site_variant', 'splice_site_variant'),
     ('genome_support', 'genome_support'),
     ]
)

def filter_events(events, min_support):
    failed = set()
    event_to_index = dict((events[i], i) for i in range(len(events)))
    for i in range(len(events)):
        if events[i].spanning < min_support:
            failed.add(i)
            if events[i].link:
                print('failed link', events[i].seq_id, events[i].event)
                if type(events[i].link) is list:
                    for e in events[i].link:
                        if e in event_to_index:
                            failed.add(event_to_index[e])
                elif events[i].link in event_to_index:
                    failed.add(event_to_index[events[i].link])
                else:
                    print('cannot find link', events[i].seq_id, events[i].event, events[i].link)

    for i in sorted(list(failed), reverse=True):
        del events[i]

def extract_features(gtfs, feature_types=('exon', 'junction')):
    annotated_features = {}
    for feature_type in feature_types:
        annotated_features[feature_type] = set()
    for gtf in gtfs:
        features = Transcript.extract_features(gtf)
        for feature_type in feature_types:
            annotated_features[feature_type] = annotated_features[feature_type].union(
                features[feature_type])
    return annotated_features


def find_novel_junctions(matches, align, transcript, query_seq, ref_fasta, accessory_known_features=None, max_diff=1):
    """Find novel junctions within a single gene/transcript

    Args:
        block_matches: (list) dictionaries where
                                          key=transcript name,
                                          value=[match1, match2, ...] where
                                                match1 = matches of each alignment block
                                                         i.e.
                                                         [(exon_id, '=='), (exon_id, '==')]
        align: (Alignment) alignment object
        transcripts_dict: (dictionary) key=transcript_id value=Transcript object
        ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
        accessary_features: (dict) 'exon' or 'junction': tuples of (chrom, start, end)
    Returns:
        List of event (dictionary storing type of event, exon indices, genomic coordinates and size)
    """
    def event_to_adj(event):
        adj = Adjacency(align.query,
                        (align.target, align.target),
                        (event['seq_breaks'][0], event['seq_breaks'][1]),
                        (event['pos'][0], event['pos'][1]),
                        transcripts=(transcript, transcript),
                        chroms=(align.target, align.target),
                        genome_breaks=(event['pos'][0], event['pos'][1]),
                        event=event['event'],
                        size=event['size'],
                        )

        if event['event'] == 'skipped_exon':
            adj.exons = (transcript.coord_to_exon(event['pos'][0]),
                         transcript.coord_to_exon(event['pos'][1]))

            if event['blocks']:
                # set up splice_motif field for genome association
                check_splice_motif_skipped_exon(adj, transcript, ref_fasta)

        elif event['event'] == 'novel_exon':
            exon1 = transcript.coord_to_exon(event['pos'][0])
            exon2 = transcript.coord_to_exon(event['pos'][1])
            if exon1 is None and exon2 is not None:
                adj.exons = ('na', exon2)
            elif exon1 is not None and exon2 is None:
                adj.exons = (exon1, 'na')

        else:
            if len(event['exons']) == 1:
                adj.exons = (transcript.exon_num(event['exons'][0]),
                             transcript.exon_num(event['exons'][0]))
            else:
                adj.exons = (transcript.exon_num(event['exons'][0]),
                             transcript.exon_num(event['exons'][1]))

        if 'splice_motif' in event:
            adj.splice_motif = event['splice_motif']

        if adj.event == 'retained_intron':
            adj.set_probe(query_seq)
        else:
            adj.set_probe(query_seq)

        set_frame(adj)

        if adj.event in ('novel_acceptor', 'novel_donor') and adj.splice_motif[1] is None:
            check_splice_motif_novel_site(adj, transcript, ref_fasta)

        return adj

    def set_frame(adj):
        if adj.size % 3 == 0:
            adj.in_frame = True
        else:
            adj.in_frame = False

    def split_event(event):
        event1 = deepcopy(event)
        event2 = deepcopy(event)

        if event['event'] == 'retained_intron':
            event1['pos'] = (event['pos'][0], event['pos'][0] + 1)
            event1['seq_breaks'] = (event['seq_breaks'][0], event['seq_breaks'][0] + 1)
            event2['pos'] = (event['pos'][1] - 1, event['pos'][1])
            event2['seq_breaks'] = (event['seq_breaks'][1] - 1, event['seq_breaks'][1])
        else:
            event1['pos'] = (align.blocks[event['blocks'][0] - 1][1], event['pos'][0])
            event1['seq_breaks'] = (event['seq_breaks'][0] - 1, event['seq_breaks'][0])
            event2['pos'] = (event['pos'][1], align.blocks[event['blocks'][-1] + 1][0])
            event2['seq_breaks'] = (event['seq_breaks'][1], event['seq_breaks'][1] + 1)

        adj1 = event_to_adj(event1)
        adj2 = event_to_adj(event2)
        adj1.link = adj2
        adj2.link = adj1
        return (adj1, adj2)

    def remove_redundant_splice_sites_of_novel_intron(events):
        """Remove novel_acceptor and novel_donor corresponding to novel_intron"""
        novel_introns = [e for e in events if e['event'] == 'novel_intron']

        for ni in novel_introns:
            acceptor_index = None
            donor_index = None
            for i in range(len(events)):
                if ni['pos'] == events[i]['pos']:
                    if events[i]['event'] == 'novel_acceptor':
                        acceptor_index = i
                    elif events[i]['event'] == 'novel_donor':
                        donor_index = i
            if acceptor_index is not None and donor_index is not None:
                for i in sorted([acceptor_index, donor_index], reverse=True):
                    del events[i]

    # find annotated junctions
    annotated = set()
    # sort multiple exon matches for single exon by exon num
    [m.sort(key=lambda mm: int(mm[0])) for m in matches if m is not None]
    for i in range(len(matches) - 1):
        j = i + 1

        if matches[i] is None or matches[j] is None:
            continue

        if (i, j) in annotated:
            continue

        if is_junction_annotated(matches[i][-1], matches[j][0]):
            annotated.add((i, j))
            continue

    known_exons = known_juncs = None
    if accessory_known_features is not None:
        known_exons = accessory_known_features['exon']
        known_juncs = accessory_known_features['junction']

    all_events = []
    j = 0
    for i in range(len(matches)):
        # for retained intron, not a 'junction'
        if matches[i] is not None and len(matches[i]) > 1:
            events = classify_novel_junction(matches[i],
                                             None,
                                             align.target,
                                             align.blocks[i],
                                             transcript,
                                             ref_fasta,
                                             known_juncs=known_juncs,
                                             known_exons=known_exons,
                                             )
            if events:
                for e in events:
                    e['blocks'] = (i, j)
                    e['transcript'] = transcript.id
                    e['seq_breaks'] = [align.query_blocks[i][0], align.query_blocks[i][1]]
                all_events.extend(events)

        if i == len(matches) - 1:
            break

        j = i + 1

        if matches[i] is None and matches[j] is not None:
            # special case where the first 2 blocks is the utr and there's an
            # insertion separating the 2 blocks
            if i == 0:
                events = classify_novel_junction(matches[i],
                                                 matches[j][0],
                                                 align.target,
                                                 align.blocks[i:j + 1],
                                                 transcript,
                                                 ref_fasta,
                                                 known_juncs=known_juncs,
                                                 known_exons=known_exons,
                                                 )
                for e in events:
                    e['blocks'] = (i, j)
                    e['transcript'] = transcript.id
                    e['seq_breaks'] = [align.query_blocks[i][1], align.query_blocks[j][0]]
                all_events.extend(events)
            continue

        # skip if junction is annotated
        if (i, j) in annotated:
            continue

        # for novel exon
        if matches[j] is None and matches[i] is not None:
            for k in range(j + 1, len(matches)):
                if matches[k] is not None and matches[k][0][0] - matches[i][-1][0] == 1:
                    j = k
                    break

        if matches[i] is not None and matches[j] is not None:
            # matches (i and j) are the flanking matches, blocks are the middle
            # "novel" blocks
            events = classify_novel_junction(matches[i][-1],
                                             matches[j][0],
                                             align.target,
                                             align.blocks[i:j + 1],
                                             transcript,
                                             ref_fasta,
                                             known_juncs=known_juncs,
                                             known_exons=known_exons,
                                             max_diff=max_diff,
                                             )

            if events:
                for e in events:
                    e['blocks'] = range(i + 1, j)
                    e['seq_breaks'] = [align.query_blocks[i][1], align.query_blocks[j][0]]
                    e['transcript'] = transcript.id
                all_events.extend(events)

    # novel acceptor and donor corresponding to novel_intron are called, remove
    remove_redundant_splice_sites_of_novel_intron(all_events)

    adjs = []
    for event in all_events:
        if event['event'] in ('novel_exon', 'retained_intron'):
            if event['event'] == 'retained_intron':
                if align.strand == '+':
                    event['seq_breaks'][0] += event['seq_break_offsets'][0]
                    event['seq_breaks'][1] -= event['seq_break_offsets'][1]
                else:
                    event['seq_breaks'][0] -= event['seq_break_offsets'][0]
                    event['seq_breaks'][1] += event['seq_break_offsets'][1]

            adjs.extend(split_event(event))

            if event['event'] == 'retained_intron':
                check_splice_motif_retained_intron(event, adjs[-2:], query_seq, ref_fasta, align, transcript, max_diff=max_diff)
        else:
            adjs.append(event_to_adj(event))

    return adjs


def report(event, event_id=None):
    data = []

    data.append(event.chroms[0])
    data.append(event.genome_breaks[0] - 1)
    data.append(event.genome_breaks[0])
    data.append(event.chroms[0])
    data.append(event.genome_breaks[1])
    data.append(event.genome_breaks[1] + 1)
    data.append(event_id)
    data.append('.')
    data.append('+')
    data.append('+')

    for item, label in report_items.items():
        value = 'na'

        if label is not None and ',' in label:
            attr, index = label.split(',')
            if hasattr(event, attr):
                values = getattr(event, attr)
                if (isinstance(values, tuple) or isinstance(values, list)) and\
                   len(values) > int(index):
                    value = values[int(index)]

        elif item == 'gene':
            value = event.transcripts[0].gene

        elif item == 'transcript':
            value = event.transcripts[0].id

        elif item == 'splice_site_variant' and value == 'na':
            if hasattr(event, 'splice_site_variant'):
                value = getattr(event, 'splice_site_variant')

            elif hasattr(event, 'splice_motif'):
                val = getattr(event, 'splice_motif')
                if val is not None and val[0] is not None and len(val) > 0 and val[1] is not None and len(val[1][0]) != 2:
                    variants = []
                    for v in val[1]:
                        variants.append(v[2])
                    value = ','.join(variants)

        elif hasattr(event, label):
            val = getattr(event, label)

            if label == 'splice_motif' and val is not None:
                val = val[0]

            if val is not None:
                value = val

        data.append(str(value))

    return '\t'.join(map(str, data))


def classify_novel_junction(match1, match2, chrom, blocks, transcript, ref_fasta, min_intron_size=20,
                            known_juncs=None, known_exons=None, no_indel=True, max_diff=1):
    """Classify given junction into different splicing events or indel

    Args:
        match1: (tuple or list) single tuple: exon_index, 2-char match result e.g. '==', '>=', etc
                                list of 2 tuples for retained_intron (special case)
        match2: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
        chrom: (str) chromosome, for checking splice motif
        blocks: (list) list of list (block coordinates)
        transcript: (Transcript) object of transcript, for getting exon coordinates
        ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
    Returns:
        List of event (dictionary storing type of event, exon indices, genomic coordinates and size)
    """
    events = []
    # set default values for 'pos'
    if isinstance(blocks[0], int):
        pos = (blocks[0], blocks[1])
        if known_exons and (chrom, pos[0], pos[1]) in known_exons:
            return events
    else:
        pos = (blocks[0][1], blocks[1][0])
        if known_juncs and (chrom, pos[0], pos[1]) in known_juncs:
            return events

    if match2 is None:
        if len(match1) == 2:
            exons = [m[0] for m in match1]
            # only requires one flanking exon boundary matched, on the other
            # side it's OK if the contig does not reach
            if ((match1[0][1] == '=>' and match1[-1][1][0] == '<') or
                (match1[0][1][1] == '>' and match1[-1][1] == '<=')) and\
               len([(a, b) for a, b in zip(exons, exons[1:]) if b == a + 1]) == len(match1) - 1:
                size = transcript.exons[exons[1]][0] - transcript.exons[exons[0]][1] - 1
                if not known_exons or not (chrom, pos[0], pos[1]) in known_exons:
                    events.append({'event': 'retained_intron',
                                   'exons': exons,
                                   'pos': (transcript.exons[exons[0]][1], transcript.exons[exons[1]][0]),
                                   'seq_break_offsets': (transcript.exons[exons[0]][1] - pos[0],
                                                         pos[1] - transcript.exons[exons[1]][0]),
                                   'size': size})

        # a single block spans multiple exons e.g. [(4, '=>'), (5, '<>'), (6,
        # '<<')]
        else:
            retained_introns = []
            new_event = None
            for i in range(len(match1) - 1):
                m1 = match1[i]
                m2 = match1[i + 1]
                if m2[0] == m1[0] + 1 and m1[1][1] == '>' and m2[1][0] == '<':
                    exons = m1[0], m2[0]
                    size = transcript.exons[exons[1]][0] - transcript.exons[exons[0]][1] - 1
                    if not known_exons or not (chrom, pos[0], pos[1]) in known_exons:
                        events.append({'event': 'retained_intron',
                                       'exons': exons,
                                       'pos': (transcript.exons[exons[0]][1], transcript.exons[exons[1]][0]),
                                       'seq_break_offsets': (transcript.exons[exons[0]][1] - pos[0],
                                                             pos[1] - transcript.exons[exons[1]][0]),
                                       'size': size})

        # check if retained_intron is known
        already_known_indices = []
        if known_exons:
            for i in range(len(events)):
                if events[i]['event'] == 'retained_intron':
                    encompassing_exons = [e for e in known_exons if e[0] == chrom and
                                          e[1] <= events[i]['pos'][0] and
                                          e[2] >= events[i]['pos'][1]]
                    if encompassing_exons:
                        already_known_indices.append(i)
            for i in reversed(already_known_indices):
                del events[i]

    # genomic deletion or novel_intron
    elif match1 is None and match2 is not None:
        if match2[0] == 0 and match2[1] == '<=':
            if transcript.strand == '+':
                donor_start = blocks[0][1] + 1
                acceptor_start = blocks[1][0] - 2
            else:
                donor_start = blocks[1][0] - 2
                acceptor_start = blocks[0][1] + 1

            gap_size = blocks[1][0] - blocks[0][1] - 1
            pos = (blocks[0][1], blocks[1][0])
            event = None
            if gap_size > 0:
                splice_motif = check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta, max_diff=max_diff)
                if gap_size > min_intron_size and splice_motif[0]:
                    event = 'novel_intron'
                else:
                    event = 'del'
            if event is not None:
                events.append({'event': event,
                               'exons': [match2[0]],
                               'pos': pos,
                               'size': gap_size,
                               'splice_motif': splice_motif})

    else:
        if match2[0] > match1[0] + 1:
            # if match2[0] > match1[0] + 1 and\
           # '=' in match1[1] and\
           # '=' in match2[1]:
            size = 0
            for e in range(match1[0] + 1, match2[0]):
                exon = transcript.exons[e]
                size += exon[1] - exon[0] + 1
            events.append({'event': 'skipped_exon', 'exons': range(match1[0] + 1, match2[0]), 'pos': pos, 'size': size})

        if match1[0] == match2[0] and match1[1][1] == '<' and match2[1][0] == '>':
            if transcript.strand == '+':
                donor_start = blocks[0][1] + 1
                acceptor_start = blocks[1][0] - 2
            else:
                donor_start = blocks[1][0] - 2
                acceptor_start = blocks[0][1] + 1

            gap_size = blocks[1][0] - blocks[0][1] - 1
            pos = (blocks[0][1], blocks[1][0])
            event = None
            splice_motif = [None]
            if gap_size > 0:
                splice_motif = check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta, max_diff=max_diff)
                if gap_size > min_intron_size and splice_motif[0]:
                    event = 'novel_intron'
                else:
                    event = 'del'
            elif gap_size == 0:
                event = 'ins'

            if event is not None:
                if event != 'ins':
                    events.append({'event': event,
                                   'exons': [match1[0]],
                                   'pos': pos,
                                   'size': gap_size,
                                   'splice_motif': splice_motif})
                else:
                    events.append({'event': event, 'exons': [match1[0]], 'pos': pos, 'splice_motif': splice_motif})

        # novel donor and acceptor
        if match2[1][0] != '=':
            # if match1[1][1] == '=' and match2[1][0] != '=':
            size = abs(blocks[1][0] - transcript.exons[match2[0]][0])
            if transcript.strand == '+':
                event = 'novel_acceptor'
                donor_start = blocks[0][1] + 1
                acceptor_start = blocks[1][0] - 2
            else:
                event = 'novel_donor'
                donor_start = blocks[1][0] - 2
                acceptor_start = blocks[0][1] + 1
            # check splice motif
            splice_motif = check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta, max_diff=max_diff)
            if splice_motif[0]:
                events.append({'event': event,
                               'exons': [match1[0], match2[0]],
                               'pos': pos,
                               'size': size,
                               'splice_motif': splice_motif})

        if match1[1][1] != '=':
            # if match1[1][1] != '=' and match2[1][0] == '=':
            size = abs(blocks[0][1] - transcript.exons[match1[0]][1])
            if transcript.strand == '+':
                event = 'novel_donor'
                donor_start = blocks[0][1] + 1
                acceptor_start = blocks[1][0] - 2
            else:
                event = 'novel_acceptor'
                donor_start = blocks[1][0] - 2
                acceptor_start = blocks[0][1] + 1
            # check splice motif
            splice_motif = check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta, max_diff=max_diff)
            if splice_motif[0]:
                events.append({'event': event,
                               'exons': [match1[0], match2[0]],
                               'pos': pos,
                               'size': size,
                               'splice_motif': splice_motif})

        if match2[0] == match1[0] + 1 and\
           match1[1][1] == '=' and\
           match2[1][0] == '=':
            pos = (blocks[1][0], blocks[-2][1])
            if not known_exons or not (chrom, pos[0], pos[1]) in known_exons:
                size = blocks[-2][1] - blocks[1][0] + 1
                if transcript.strand == '+':
                    donor_start = pos[1] + 1
                    acceptor_start = pos[0] - 2
                else:
                    donor_start = pos[0] - 2
                    acceptor_start = pos[1] + 1

                # check splice motif
                splice_motif = check_splice_motif(chrom, donor_start, acceptor_start, transcript.strand, ref_fasta, max_diff=max_diff)
                if splice_motif[0]:
                    events.append({'event': 'novel_exon',
                                   'exons': [],
                                   'pos': pos,
                                   'size': size,
                                   'splice_motif': splice_motif})

    # set size to None for event that doesn't have size i.e. 'ins'
    for event in events:
        if 'size' not in event:
            event['size'] = None

    if no_indel:
        return [e for e in events if e['event'] not in ('ins', 'del')]
    else:
        return events


def find_diff(seq1, seq2):
    diff = []
    i = 0
    for b1, b2 in zip(seq1.lower(), seq2.lower()):
        if b1 != b2:
            diff.append((i, b1, b2))
        i += 1

    return diff


def check_splice_motif(chrom, donor_start, acceptor_start, strand, ref_fasta, max_diff=1):
    """Curates splice motif sequence and deviation from canoncial sequence

       Args:
            chrom: (str) chromosome
            donor_start: (int) in relation to reference, not transcript
            acceptor_start: (int) in relation to reference, not transcript
            strand: (str) transcript strand, '+' or '-'
            ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence, for checking splice motif
            max_diff: maximum number of bases deviated from canonical motif
       Returns:
            2-member-tuple: (str) motif sequence observed
                            (list of tuples): list of deviations from canoncial motif
                                              each deviation = (genomic coordinate,
                                                               observed base,
                                                               variant = chrom:coordBaseFrom>BaseTo)
    """
    result = [None]

    # must be at least 1 bp separating the donor and acceptor
    if abs(acceptor_start - donor_start) < 3:
        return result

    donor_seq = ref_fasta.fetch(chrom, donor_start - 1, donor_start - 1 + 2)
    acceptor_seq = ref_fasta.fetch(chrom, acceptor_start - 1, acceptor_start - 1 + 2)

    coords = [donor_start, acceptor_start]
    seqs = [donor_seq, acceptor_seq]
    canonicals = ['gt', 'ag']
    diffs = []
    for i, seq in enumerate(seqs):
        if strand == '-':
            diff = find_diff(canonicals[i], reverse_complement(seq))
        else:
            diff = find_diff(canonicals[i], seq)
        diffs.append(diff)

    base_diffs = []
    # don't allow both bases to be mutated in donor or acceptor
    if len(diffs[0]) + len(diffs[1]) <= max_diff and\
       len(diffs[0]) <= 1 and len(diffs[1]) <= 1:
        if strand == '-':
            motif = (reverse_complement(donor_seq) + reverse_complement(acceptor_seq)).upper()
        else:
            motif = (donor_seq + acceptor_seq).upper()

        # perfect match
        if not diffs[0] and not diffs[1]:
            return [motif, None]

        # 0 = donor, 1 = acceptor
        for i, diff in enumerate(diffs):
            if not diff:
                continue
            if strand == '+':
                coord_diff = coords[i] + diff[0][0]
                base_from = diff[0][2].upper()
                base_diff = diff[0][1].upper()
                variant = '{}:{}{}>{}'.format(chrom, coord_diff, base_from, base_diff)
            else:
                coord_diff = coords[i] + 1 - diff[0][0]
                base_from = reverse_complement(diff[0][2].upper())
                base_diff = reverse_complement(diff[0][1].upper())
                variant = '{}:{}{}>{}'.format(chrom, coord_diff, base_from, base_diff)
            base_diffs.append([coord_diff, base_diff, variant])

        result = [motif, base_diffs]

    return result


def check_splice_motif_retained_intron(event, adjs, query_seq, ref_fasta, align, transcript, max_diff=1):
    """Looks for base change in splice motif captured by contig in retained_introns

       Args:
           event: (dict) retained_intron from find_novel_junctions() before splitting event into adjacencies
           adjs: (list of Adjacencies) Adjacencies for setting the splice_motif attribute
           ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence,
                      for extracting splice motif
           align: (Alignment) for extracting alignment strand in generating retained_intro sequence
                              based on genome positive strand
           transcript: (Transcript) for orienting splice motif sequences
           max_diff: (int) maximum number of base subsitutions to be found in splice_motif

       Returns: None
                It will reset splice motif of each of the 2 adjacencies to tuple:
                first: 2-base splice motif observed corresponding to each adjacency
                second: None if canoncial splice motif observed
                        list of (genoome coordinate, new base, variant) if non-canonical splice motif observed
    """
    # seq_breaks = contig breakpoints, sort for extracting contig sequence in
    # case not sorted
    seq_breaks_sorted = sorted(event['seq_breaks'])
    retained_seq_contig = query_seq[seq_breaks_sorted[0]
        :seq_breaks_sorted[1] - 1]
    if align.strand == '-':
        retained_seq_contig = reverse_complement(retained_seq_contig)
    # assume 'pos' is always sorted
    retained_seq_genome = ref_fasta.fetch(adjs[0].chroms[0], event['pos'][0], event['pos'][1] - 1)
    # splice motif from assembled sequence (in contig orientation)
    splice_motif = retained_seq_contig[:2].upper() + retained_seq_contig[-2:].upper()

    # motifs[0] = upstream, motifs[1] = downstream, regardless of
    # donor/acceptor
    motifs = []
    if transcript.strand == '+':
        motifs.append(splice_motif[:2])
        motifs.append(splice_motif[-2:])
    else:
        splice_motif = reverse_complement(splice_motif)
        motifs.append(splice_motif[-2:])
        motifs.append(splice_motif[:2])

    coords = [event['pos'][0] + 1, event['pos'][1] - 2]
    seqs = [retained_seq_contig[:2], retained_seq_contig[-2:]]
    canonicals = [retained_seq_genome[:2], retained_seq_genome[-2:]]
    diffs = []
    for i, seq in enumerate(seqs):
        diff = find_diff(canonicals[i], seq)
        diffs.append(diff)

    if len(diffs[0]) + len(diffs[1]) <= max_diff and\
       len(diffs[0]) <= 1 and len(diffs[1]) <= 1:
        # 0 = donor, 1 = acceptor
        for i, diff in enumerate(diffs):
            if not diff:
                adjs[i].splice_motif = motifs[i], None
                continue

            coord_diff = coords[i] + diff[0][0]
            base_from = diff[0][1].upper()
            base_diff = diff[0][2].upper()
            variant = '{}:{}{}>{}'.format(adjs[0].chroms[0], coord_diff, base_from, base_diff)
            adjs[i].splice_motif = motifs[i], [
                [coord_diff, base_diff, variant]]


def check_splice_motif_skipped_exon(adj, transcript, ref_fasta):
    """Gathers reference exon splice site coordinates and bases for searching for
       canonical splice-site disruption - skipped exons

       only checks for splice-motifs in introns

       Args:
           adj: (Adjacency) skipped_exon
           transcript: (Transcript) transcript object for extracting strand, chromosome
           ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence,
                      for extracting splice motif
       Returns: None
               It will reset splice_motif attribute to a list of 2
               first: splice_motif sequence (should be canonical), will arrange bases according
                      to transcript strand
               second: list of [coordinate, base]
                       (each member of adj.splice_motif[1] should have 2 members)
    """
    # identify flanking splice-site coordinates of skipped exon(s)
    flanking_exons = sorted(adj.exons)
    skipped_exons = sorted(
        set(xrange(min(flanking_exons), max(flanking_exons))).difference(flanking_exons))

    exon_coords = set()
    for exon in skipped_exons:
        bounds = transcript.exon(exon)
        exon_coords.add(bounds[0])
        exon_coords.add(bounds[1])
    sorted_exon_coords = sorted(exon_coords)
    splice_site_coords = range(
        sorted_exon_coords[0] - 2,
        sorted_exon_coords[0]) + range(
        sorted_exon_coords[1] + 1,
        sorted_exon_coords[1] + 3)

    # extract coordinates and corresponding bases for
    # find_genome_interruptions()
    adj.splice_motif = [None, []]
    for coord in splice_site_coords:
        adj.splice_motif[1].append([coord, ref_fasta.fetch(transcript.chrom, coord - 1, coord).upper()])

    # fills in splice_motif sequence
    adj.splice_motif[0] = ''.join([b[1] for b in adj.splice_motif[1]])
    if transcript.strand == '-':
        adj.splice_motif[0] = reverse_complement(adj.splice_motif[0][-2:] + adj.splice_motif[0][:2])

def check_splice_motif_novel_site(adj, transcript, ref_fasta):
    """Gathers reference exon splice site coordinates and bases for searching for
       canonical splice-site disruption - novel_donor/acceptor

       only checks for splice-motifs in introns

       Args:
           adj: (Adjacency) novel_donor/novel_acceptor
                second member of splice_motif attribute (list) should be None
           transcript: (Transcript) transcript object for extracting chromosome
           ref_fasta: (Pysam.Fastafile) Pysam handle to access reference sequence,
                      for extracting splice motif
      Returns: None
               Sets the second member of splice_motif attribute to list of [coordinate, base]
               (each member of adj.splice_motif[1] should have 2 members)
    """
    splice_site_coords = []
    for exon in adj.exons:
        splice_site_coords.extend(transcript.exon(exon))
    sorted_splice_site_coords = sorted(splice_site_coords)

    adj.splice_motif[1] = []
    for i in range(1, 3):
        if i == 1:
            for j in range(1, 3):
                coord = sorted_splice_site_coords[i] + j
                adj.splice_motif[1].append([coord, ref_fasta.fetch(transcript.chrom, coord - 1, coord).upper()])
        if i == 2:
            for j in range(2, 0, -1):
                coord = sorted_splice_site_coords[i] - j
                adj.splice_motif[1].append([coord, ref_fasta.fetch(transcript.chrom, coord - 1, coord).upper()])


def is_junction_annotated(match1, match2):
    """Checks if junction is in gene model

    Args:
        match1: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
        match2: (tuple) exon_index, 2-char match result e.g. '==', '>=', etc
    Returns:
        True or False
    """
    if match2[0] == match1[0] + 1 and match1[1][1] == '=' and match2[1][0] == '=':
        return True

    return False


def corroborate_genome(adjs, genome_bam=False, vcfs=[], sample=None):
    """Check if genome has splice-site mutation leading to novel splice site

    Args:
        adjs: list of adjs
        genome_bam: genome bam
    """
    canonical = 'GTAG'
    for adj in adjs:
        if adj.splice_motif is None:
            continue

        if adj.event in ('novel_acceptor', 'novel_donor', 'retained_intron', 'novel_intron', 'novel_exon'):
            if adj.splice_motif[0] != canonical and adj.splice_motif[1]:
                if genome_bam:
                    find_genome_support(adj, genome_bam)
                elif vcfs:
                    find_vcf_support(adj, vcfs, sample=sample)

        if adj.event in ('skipped_exon', 'novel_acceptor', 'novel_donor'):
            if adj.splice_motif[0] == canonical and adj.splice_motif[1] and len(adj.splice_motif[1]) == 4:
                if genome_bam:
                    find_genome_interruptions(adj, genome_bam)
                elif vcfs:
                    find_vcf_interruptions(adj, vcfs, sample=sample)

def find_genome_support(adj, bam):
    """Reports level of support of novel base change from genome bam file

    Args:
        adj: (Adjacency) novel_acceptor/donor (non-canonical splice motif in reference observed)
                         retained_intron (non-canonical splice motif captured in contig)
        bam: (Pysam bam handle) genome bam
    """
    chrom = adj.chroms[0]
    if not adj.chroms[0] in bam.references:
        if adj.chroms[0][:3] == 'chr':
            chrom = adj.chroms[0].lstrip('chr')
        else:
            chrom = 'chr' + adj.chroms[0]
        if chrom not in bam.references:
            return None

    genome_support = []
    for pos, base, variant in adj.splice_motif[1]:
        counts = {'corroborated': 0, 'total': 0}
        for pileup_col in bam.pileup(chrom, pos - 2, pos - 1):
            for pileup_read in pileup_col.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    if pileup_col.pos == pos - 1:
                        counts['total'] += 1
                        if pileup_read.alignment.query_sequence[pileup_read.query_position].upper() == base.upper():
                            counts['corroborated'] += 1

        genome_support.append('{}/{}'.format(counts['corroborated'], counts['total']))

    adj.genome_support = ','.join(genome_support)

def extract_support_vcf_rec(rec, sample=None):
    ''' assumptions: first allele is ref, second allele is alt, only one alt allele, DP only single number '''
    if sample is None:
        DPs = [s['DP'] for s in rec.samples.values()]
        ADs = [s['AD'] for s in rec.samples.values()]
        return '{}/{}'.format(ADs[0][1], DPs[0])
    else:
        for s in rec.samples:
            if s == sample:
                return '{}/{}'.format(rec.samples[s]['AD'][1], rec.samples[s]['DP'])

def find_vcf_support(adj, vcfs, sample=None):
    genome_support = []
    for pos, base, variant in adj.splice_motif[1]:
        for vcf in vcfs:
            for rec in vcf.fetch(adj.chroms[0], pos-1, pos):
                if base.upper() in rec.alleles and pos == rec.pos:
                    support = extract_support_vcf_rec(rec, sample)
                    if support:
                        genome_support.append(support)

    if genome_support:
        adj.genome_support = ','.join(genome_support)

def find_genome_interruptions(adj, bam):
    """Search for genomic variants that disrupt splice motifs

    Args:
        adj: (Adjacency) novel_acceptor/donor (reference splice junction not used, novel junction canonical)
                         skipped_exon
        bam: (Pysam bam handle) genome bam
    """
    chrom = adj.chroms[0]
    if not adj.chroms[0] in bam.references:
        if adj.chroms[0][:3] == 'chr':
            chrom = adj.chroms[0].lstrip('chr')
        else:
            chrom = 'chr' + adj.chroms[0]
        if chrom not in bam.references:
            return None

    variants = []
    supports = []
    for pos, ref_base in adj.splice_motif[1]:
        counts = {}
        for pileup_col in bam.pileup(chrom, pos - 2, pos - 1):
            for pileup_read in pileup_col.pileups:
                if not pileup_read.is_del and not pileup_read.is_refskip:
                    if pileup_col.pos == pos - 1:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        if base not in counts:
                            counts[base] = 0
                        counts[base] += 1

        sorted_bases = sorted(counts.keys(), key=lambda b: counts[b], reverse=True)
        coverage = sum(counts.values())

        for base in sorted_bases:
            if base.upper() == ref_base.upper():
                continue
            if counts[base] > coverage / 4:
                variants.append('{}:{}{}>{}'.format(adj.chroms[0], pos, ref_base.upper(), base.upper()))
                supports.append('{}/{}'.format(counts[base], coverage))
            break

    if variants:
        adj.splice_site_variant = ','.join(variants)
        adj.genome_support = ','.join(supports)

def find_vcf_interruptions(adj, vcfs, sample=None, search_dist=500):
    variants = []
    supports = []
    for pos, ref_base in adj.splice_motif[1]:
        for vcf in vcfs:
            for rec in vcf.fetch(adj.chroms[0], pos-search_dist, pos+search_dist):
                if rec.pos == pos and rec.rlen == 1:
                    variants.append('{}:{}{}>{}'.format(adj.chroms[0], pos, ref_base.upper(), ','.join([a.upper() for a in list(rec.alts)])))
                    support = extract_support_vcf_rec(rec, sample)
                    if support:
                        supports.append(support)
                elif rec.start <= pos and rec.stop >= pos:
                    adj.splice_site_variant = ','.join(list(rec.alts))
                    variants.append('{}:{}{}>{}'.format(adj.chroms[0], '{}-{}'.format(rec.start, rec.stop), rec.ref, ','.join([a.upper() for a in list(rec.alts)])))
                    support = extract_support_vcf_rec(rec, sample)
                    if support:
                        supports.append(support)

    if variants:
        adj.splice_site_variant = ','.join(variants)
        if supports:
            adj.genome_support = ','.join(supports)
