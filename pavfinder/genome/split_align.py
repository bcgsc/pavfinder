import pysam
from intspan import intspan
import sys
import re
from . import bwa_mem
from .adjacency import Adjacency
from .alignment import compare_chr, Alignment, target_non_canonical


def find_chimera(alns, bam, min_coverage=0.95, check_alt_paths=False, max_splits=3, check_haplotype=True, debug=False):
    """Finds primary_aligns alignments corresponding to chimera

    Args:
        bam: Pysam handle to bam file
    Returns:
        List of Alignments
    """

    def _is_secondary_chimera_better(primary_chimera, secondary_chimera):
        primary_coverage = get_contig_coverage(primary_chimera)
        secondary_coverage = get_contig_coverage(secondary_chimera)

        if secondary_coverage > primary_coverage:
            return True
        elif secondary_coverage == primary_coverage:
            primary_edit_distance = 0
            secondary_edit_distance = 0
            for align in primary_chimera:
                primary_edit_distance += bwa_mem.effective_edit_distance(align.sam)
            for align in secondary_chimera:
                secondary_edit_distance += bwa_mem.effective_edit_distance(align.sam)
            if secondary_edit_distance <= primary_edit_distance:
                return True
            else:
                return False
        else:
            return False

    def _has_end_to_end_secondary_align(aligns):
        for align in aligns:
            ops = [a[0] for a in align.sam.cigar]
            if ops[0] == 0 and ops[-1] == 0:
                return align
        return None

    def _replace_trl(primary_aligns, secondary_aligns, max_span=1000000):
        alt_paths = []
        for primary_align in primary_aligns:
            aligns = [align for align in secondary_aligns if align.target == primary_align.target] + [primary_align]
            alt_paths = [path for path in find_paths(aligns, min_coverage=min_coverage, max_ends=50, get_all=True, debug=debug)
                         if len(path) == 2 and abs(aligns[path[0]].tstart - aligns[path[1]].tstart) < max_span]

            if alt_paths:
                alt_paths.sort(key=lambda path: abs(aligns[path[0]].tstart - aligns[path[1]].tstart))
                return [aligns[i] for i in alt_paths[0]]

        return None

    primary_aligns, secondary_aligns = bwa_mem.find_chimera(alns, bam, check_haplotype=check_haplotype, debug=debug)

    if primary_aligns:
        if len(primary_aligns) > max_splits:
            if debug:
                sys.stdout.write('{} number of primary split alignments({}) exceeds maximum:{}\n'.format(primary_aligns[0].query,
                                                                                                         len(primary_aligns),
                                                                                                         max_splits))
            return None, None

        primary_paths = find_paths(primary_aligns, min_coverage=min_coverage, debug=debug)

        if primary_paths:
            if secondary_aligns:
                align = _has_end_to_end_secondary_align(secondary_aligns)
                if align:
                    if debug:
                        sys.stdout.write('filter out chimera {} ({}bp) because of end-to-end secondary alignment {}:{}-{} {}\n'.format(align.query,
                                                                                                                                       align.query_len,
                                                                                                                                       align.target,
                                                                                                                                       align.tstart,
                                                                                                                                       align.tend,
                                                                                                                                       align.sam.cigarstring))
                    return None, None

            primary_chimera = [primary_aligns[i] for i in primary_paths]
            primary_chroms = set([align.target for align in primary_chimera])
            for align in primary_chimera:
                align.dubious = False

            primary_replaced = False
            if len(primary_chroms) > 1 and secondary_aligns:
                replaced_chimera = _replace_trl(
                    primary_aligns, secondary_aligns)
                if replaced_chimera:
                    if debug:
                        sys.stdout.write('replaced {} trl {}:{}-{} {}:{}-{} with {}:{}-{} {}:{}-{}\n'.format(alns[0].qname,
                                                                                                             primary_chimera[0].target,
                                                                                                             primary_chimera[0].tstart,
                                                                                                             primary_chimera[0].tend,
                                                                                                             primary_chimera[1].target,
                                                                                                             primary_chimera[1].tstart,
                                                                                                             primary_chimera[1].tend,
                                                                                                             replaced_chimera[0].target,
                                                                                                             replaced_chimera[0].tstart,
                                                                                                             replaced_chimera[0].tend,
                                                                                                             replaced_chimera[1].target,
                                                                                                             replaced_chimera[1].tstart,
                                                                                                             replaced_chimera[1].tend,))
                    primary_replaced = True
                    primary_chimera = replaced_chimera
                    primary_chroms = set([align.target for align in primary_chimera])

            if debug:
                for align in primary_chimera:
                    output_line = ' '.join(map(str, [align.query,
                                                     align.qstart,
                                                     align.qend,
                                                     align.target,
                                                     align.tstart,
                                                     align.tend,
                                                     align.sam.cigarstring,
                                                     align.strand]))
                    sys.stdout.write('{}\n'.format(output_line))

            # remove alignments to non-canonical chromosomes
            for align in secondary_aligns[:]:
                if target_non_canonical(align.target):
                    secondary_aligns.remove(align)

            secondary_paths = find_paths(secondary_aligns, min_coverage=min_coverage, get_all=True, debug=debug)

            # alternative chimera combining primary and secondary aligns: won't
            # work for 3 ways
            dubious = set()

            if check_alt_paths:
                for i in range(len(primary_chimera)):
                    aligns = []
                    for j in range(len(primary_chimera)):
                        if j != i:
                            aligns.append(primary_chimera[j])
                    aligns.extend(secondary_aligns)
                    alt_paths_candidates = find_paths(aligns, min_coverage=min_coverage, get_all=True, no_trim=[0], debug=debug)

                    if alt_paths_candidates:
                        if debug:
                            sys.stdout.write('{} primary alignment {}:{}-{} is dubious, can be replaced by {} secondary alignments\n'.format(primary_chimera[i].query,
                                                                                                                                             primary_chimera[i].target,
                                                                                                                                             primary_chimera[i].tstart,
                                                                                                                                             primary_chimera[i].tend,
                                                                                                                                             len(alt_paths_candidates)))
                        dubious.add(i)

            # use secondary alignments to filter primary split alignments
            if secondary_paths:
                passed = True
                for indices in secondary_paths:
                    secondary_chimera = [secondary_aligns[i] for i in indices]
                    if _is_secondary_chimera_better(
                            primary_chimera, secondary_chimera):
                        if debug:
                            sys.stdout.write('Filter out {} as alternative chimera can be formed in seconday alignments\n'.format(alns[0].qname))
                        passed = False
                        break

                if passed:
                    return primary_chimera, dubious

            else:
                return primary_chimera, dubious

    return None, None


def screen_subseq_alns(adj_aligns, subseq_alns, realign_bam, name_sep, debug=False):
    """Screens subsequence alignments against original primary_aligns alignments to weed out false-positives

    Args:
        adj_aligns: (list) Initial Alignment objects
        subseq_aligns: (list) AlignedRead objects of subsequence alignments of same adj
        realign_bam: Pysam handle to realignment BAM file
        name_sep: (str) Character used to combined various info into query name

    Returns:
        True = keep, False = discard
    """
    ambiguous_NM = 5
    passed = [True, True]
    if len(subseq_alns) >= len(adj_aligns):
        matched = {0: None, 1: None}
        alt_alns = {0: [], 1: []}
        for aln in subseq_alns:
            # check this, otherwise getrname() will fail
            # if alignment target is not canonical
            if aln.is_unmapped or re.search('[-_.]', realign_bam.getrname(aln.tid)):
                continue

            idx = int(aln.qname.split(name_sep)[-1])

            # if end-to-end match
            if re.match(r'\d+.*M$', aln.cigarstring) and not re.search('[SH]', aln.cigarstring):
                if realign_bam.getrname(aln.tid) == adj_aligns[idx].target and aln.pos + 1 == adj_aligns[idx].tstart:
                    matched[idx] = aln
                else:
                    try:
                        alt_alns[idx].append(aln)
                    except BaseException:
                        alt_alns[idx] = [aln]

        for i in (0, 1):
            if matched[i] is not None:
                if alt_alns[i]:
                    for alt_aln in alt_alns[i]:
                        if int(alt_aln.opt('NM')) - int(matched[i].opt('NM')) <= ambiguous_NM:
                            if debug:
                                sys.stdout.write('ambigous subseq aln: {} {} {} {} NM:{}\n'.format(alt_aln.qname,
                                                                                                   realign_bam.getrname(alt_aln.tid),
                                                                                                   alt_aln.pos,
                                                                                                   alt_aln.cigarstring,
                                                                                                   alt_aln.opt('NM')))
                            passed[i] = False
                            break

            else:
                if debug:
                    sys.stdout.write('Failed to map subseq end-to-end {}:{}-{} {}:{}-{}\n'.format(adj_aligns[i].query,
                                                                                                  adj_aligns[i].qstart,
                                                                                                  adj_aligns[i].qend,
                                                                                                  adj_aligns[i].target,
                                                                                                  adj_aligns[i].tstart,
                                                                                                  adj_aligns[i].tend))
                passed[i] = False

    return passed


def find_paths(aligns, min_coverage=None, use_end_to_end=True, get_all=False,
               max_nodes=500, max_paths=5, max_ends=50,
               no_trim=[], same_target=None, from_edge=0.02, debug=False):
    def _find_end_points():
        starts = []
        ends = []
        from_start = max(1, int(from_edge * aligns[0].query_len))
        from_end = aligns[0].query_len - from_start + 1
        for i in range(len(aligns)):
            if aligns[i].qstart <= from_start:
                try:
                    starts.append(i)
                except BaseException:
                    starts = [i]
            if aligns[i].qend >= from_end:
                try:
                    ends.append(i)
                except BaseException:
                    ends = [i]

        for i in starts:
            if i in ends:
                if aligns[i].qstart < aligns[i].query_len - aligns[i].qend + 1:
                    ends.remove(i)
                else:
                    starts.remove(i)
                break

        return starts, ends

    def _trim_aligns(max_mappings=100):
        regions = {}
        for i in range(len(aligns)):
            if i in no_trim:
                continue

            key = '{}-{}:{}'.format(aligns[i].qstart, aligns[i].qend, aligns[i].strand)
            try:
                regions[key].append(i)
            except BaseException:
                regions[key] = [i]

        extra = []
        for key in regions.keys():
            if len(regions[key]) > max_mappings:
                for i in regions[key][max_mappings:]:
                    extra.append(i)

        if extra:
            new_aligns = []
            for i in range(len(aligns)):
                if i not in extra:
                    new_aligns.append(aligns[i])

            return new_aligns

        return None

    def _construct_graph(max_links=100, max_olap=0.2, max_gap=0.2):
        """
        Alignments must be like to be an edge
           ----> j
        ----> i
        """
        graph = {}
        for i in range(len(aligns)):
            len_i = aligns[i].qend - aligns[i].qstart + 1
            for j in range(len(aligns)):
                if i != j:
                    len_j = aligns[j].qend - aligns[j].qstart + 1
                    olap = max(0, aligns[i].qend - aligns[j].qstart + 1)
                    gap = max(0, aligns[j].qstart - aligns[i].qend - 1)
                    # only allows gap = 0, gap > 0, or partial overlap
                    if aligns[j].qend > aligns[i].qend and aligns[j].qstart > aligns[i].qstart:
                        if float(olap) / len_i >= max_olap or\
                           float(olap) / len_j >= max_olap or\
                           float(gap) / len_i >= max_gap or\
                           float(gap) / len_j >= max_gap:
                            continue
                        try:
                            graph[i].append(j)
                        except BaseException:
                            graph[i] = [j]

        return graph

    def _bfs(graph, start, end):
        """Breath-first search
        http://stackoverflow.com/questions/8922060/breadth-first-search-trace-path
        """
        paths = []
        # maintain a queue of paths
        queue = []
        # push the first path into the queue
        queue.append([start])
        while queue:
            # get the first path from the queue
            path = queue.pop(0)
            # get the last node from the path
            node = path[-1]
            # path found
            if node == end:
                # captures all paths
                paths.append(path)
                # this will return shortest path
                # return path
            # enumerate all adjacent nodes, construct a new path and push it
            # into the queue
            for adjacent in graph.get(node, []):
                new_path = list(path)
                new_path.append(adjacent)
                queue.append(new_path)

        return paths

    def _coverage(path):
        spans = [intspan('{}-{}'.format(aligns[i].qstart, aligns[i].qend)) for i in path]
        covered = spans[0]
        overlaps = []
        for i in range(1, len(spans)):
            covered = covered.union(spans[i])

            overlap = spans[i - 1].intersection(spans[i])
            if len(overlap) > 0:
                overlaps.append(overlap)

        return covered, overlaps

    def _screen(covered, overlaps):
        passed = True
        if min_coverage is not None:
            if not use_end_to_end:
                coverage = len(covered) / float(aligns[0].query_len)
            else:
                coverage = (max(covered) - min(covered) + 1) / float(aligns[0].query_len)

            if coverage < min_coverage:
                passed = False

        return passed

    def _pick_best(path_info):
        best = {'index': None, 'covered': None, 'overlapped': None}
        for i in sorted(path_info.keys()):
            covered = len(path_info[i]['covered'])
            if not path_info[i]['overlaps']:
                overlapped = 0
            else:
                overlapped = sum([len(olap) for olap in path_info[i]['overlaps']])

            if best['index'] is None or\
               covered > best['covered'] or\
               overlapped < best['overlapped'] or\
               len(path_info[i]['chroms']) < len(path_info[best['index']]['chroms']):
                best['index'] = i
                best['covered'] = covered
                best['overlapped'] = overlapped

        return best['index']

    # check promiscuity
    trimmed_aligns = _trim_aligns()
    if trimmed_aligns:
        aligns = trimmed_aligns

    # same target
    if same_target:
        aligns = [align for align in aligns if align.target == same_target]

    if not aligns:
        return []

    # construct graph
    graph = _construct_graph()

    if len(graph.keys()) > max_nodes:
        if debug:
            sys.stdout.write('{}: too many nodes({}) to construct path(max:{})\n'.format(aligns[0].query,
                                                                                         len(graph.keys()),
                                                                                         max_nodes))
        return []

    # get end points
    starts, ends = _find_end_points()
    if not starts or not ends:
        return []

    if len(starts) > max_ends:
        starts = starts[:max_ends]
    if len(ends) > max_ends:
        ends = ends[:max_ends]

    # find paths
    paths = []
    for start in starts:
        if max_paths and len(paths) >= max_paths:
            break
        for end in ends:
            paths.extend(_bfs(graph, start, end))
            if max_paths and len(paths) >= max_paths:
                break

    # screen
    path_info = {}
    for i in range(len(paths)):
        covered, overlaps = _coverage(paths[i])
        if _screen(covered, overlaps):
            chroms = set([aligns[j].target for j in paths[i]])
            path_info[i] = {'path': paths[i], 'covered': covered, 'overlaps': overlaps, 'chroms': chroms}

    if path_info:
        if get_all:
            return [paths[i] for i in path_info.keys()]
        else:
            return paths[_pick_best(path_info)]
    else:
        return []


def get_contig_coverage(aligns, end_to_end=False):
    """Coverage of the contig by the union of the primary_aligns alignments

    Args:
        aligns: (list) All Alignments constituting a chimera

    Returns:
        Fraction corresponding to coverage
    """
    span = intspan('{}-{}'.formt(aligns[0].qstart, aligns[0].qend))
    for i in range(1, len(aligns)):
        span = span.union(intspan('{}-{}'.format(aligns[i].qstart, aligns[i].qend)))

    if not end_to_end:
        return len(span) / float(aligns[0].query_len)
    else:
        return (max(span) - min(span) + 1) / float(aligns[0].query_len)


def find_adjs(aligns, contig_seq, dubious=None, debug=False):
    """Create adjs given primary_aligns alignments

    Will also determine microhomology and untemplated sequences based on information
    given by different aligners.

    Args:
        aligns: (list) primary_aligns Alignments
        contig_seq: (str) Contig sequence

    Returns:
        List of Adjacency objects
    """
    adjs = []
    for i in range(1, len(aligns)):
        homol_seq, homol_coords = bwa_mem.find_microhomology((aligns[i - 1], aligns[i]), contig_seq)
        novel_seq = bwa_mem.find_untemplated_sequence((aligns[i - 1], aligns[i]), contig_seq)
        adj = call_event(aligns[i - 1], aligns[i],
                         homol_seq=homol_seq,
                         homol_coords=homol_coords,
                         novel_seq=novel_seq,
                         contig_seq=contig_seq,
                         debug=debug)

        adj.dubious = False
        if dubious is not None and len(dubious) > 0:
            if i - 1 in dubious:
                adj.aligns[0][0].dubious = True
            if i in dubious:
                adj.aligns[0][1].dubious = True
            adj.dubious = True

        if adj is not None:
            adj.event_id = str(len(adjs) + 1)
            adjs.append(adj)

    return adjs


def call_event(align1, align2, homol_seq=None, homol_coords=None, novel_seq='-', contig_seq=None, no_sort=False, probe_side_len=25, debug=False):
    """Curates adj based on info given by primary_aligns alignments

    Args:
        align1: First Alignment object
        align2: Second Alignment object
        homol_seq: (str) Microhomology sequence in contig
        homol_coords: (tuple) Start and end coordinates of microhomology sequence in contig
        novel_seq: (str) Untemplated sequence at breakpoint
        contig_seq: (str) Contig sequence

    Returns:
        Adjacency object
    """
    # figure out breakpoints using query positions
    breaks = [None, None]
    orients = [None, None]
    contig_breaks = [None, None]
    if align1.qstart < align2.qstart:
        aligns = [align1, align2]
        breaks[0] = align1.tend if align1.strand == '+' else align1.tstart
        orients[0] = 'L' if max(align1.tstart, align1.tend) == breaks[0] else 'R'
        breaks[1] = align2.tstart if align2.strand == '+' else align2.tend
        orients[1] = 'L' if max(align2.tstart, align2.tend) == breaks[1] else 'R'
        contig_breaks = [align1.qend, align2.qstart]

    else:
        aligns = [align2, align1]
        breaks[0] = align2.tend if align2.strand == '+' else align2.tstart
        orients[0] = 'L' if max(align2.tstart, align2.tend) == breaks[0] else 'R'
        breaks[1] = align1.tstart if align1.strand == '+' else align1.tend
        orients[1] = 'L' if max(align1.tstart, align1.tend) == breaks[1] else 'R'
        contig_breaks = [align2.qend, align1.qstart]

    if not no_sort:
        if (aligns[0].target != aligns[1].target and compare_chr(aligns[0].target, aligns[1].target) > 0) or\
           (aligns[0].target == aligns[1].target and breaks[0] > breaks[1]):
            aligns.reverse()
            breaks.reverse()
            orients.reverse()

    rearrangement = None
    if aligns[0].target != aligns[1].target:
        rearrangement = 'trl'
    elif orients[0] == 'L' and orients[1] == 'L':
        rearrangement = 'inv'
    elif orients[0] == 'L' and orients[1] == 'R':
        if breaks[0] < breaks[1]:
            if breaks[0] + 1 == breaks[1]:
                # deletion of tandem duplicaton
                if contig_breaks[0] >= contig_breaks[1]:
                    rearrangement = 'del'
                    breaks = [breaks[1] + 1, breaks[0] + (contig_breaks[0] - contig_breaks[1] + 1)]
                    homol_seq = None
                    homol_coords = None
                else:
                    rearrangement = 'ins'
            else:
                # deletion with or without microhology
                rearrangement = 'del'

        elif breaks[0] > breaks[1]:
            rearrangement = 'dup'
            # if difference is breaks is small, should be tandem dup
            # if contig_breaks[0] >= contig_breaks[1]:
            #rearrangement = 'dup'
            # insertion with microhomology
            # else:
            #rearrangement = 'ins'

        # breaks[0] == breaks[1]
        else:
            if contig_breaks[0] < contig_breaks[1]:
                rearrangement = 'ins'
            else:
                # deletion of tandem duplicaton
                rearrangement = 'del'
                breaks = [breaks[1] + 1, breaks[0] + (contig_breaks[0] - contig_breaks[1] + 1)]
                homol_seq = None
                homol_coords = None

    elif orients[0] == 'R' and orients[1] == 'R':
        rearrangement = 'inv'
    elif orients[0] == 'R' and orients[1] == 'L':
        if breaks[0] == breaks[1]:
            rearrangement = 'ins'
        else:
            rearrangement = 'dup'

    adj = None
    if rearrangement is not None:
        # use tuple for 'chroms', 'breaks', 'orients': don't want them to be
        # changed (immutable)
        adj = Adjacency((aligns[0].target, aligns[1].target),
                        (breaks[0], breaks[1]),
                        contig=align1.query,
                        contig_breaks=contig_breaks,
                        contig_sizes=int(aligns[0].query_len),
                        rearrangement=rearrangement,
                        orients=(orients[0], orients[1]),
                        homol_seq=homol_seq,
                        homol_coords=homol_coords,
                        novel_seq=novel_seq,
                        aligns=aligns,
                        align_types='split',
                        )

        if contig_seq is not None:
            adj.probes.append(Adjacency.extract_probe(contig_seq, contig_breaks, len_on_each_side=probe_side_len)[0])

    elif debug:
        sys.stdout.write(
            "cannot figure out event of primary_aligns alignment contig:{} targets:{},{} orients:{} breaks:{} contig_breaks:{}\n".format(aligns[0].query,
                                                                                                                                         aligns[0].target,
                                                                                                                                         aligns[1].target,
                                                                                                                                         orients,
                                                                                                                                         breaks,
                                                                                                                                         contig_breaks))
    return adj
