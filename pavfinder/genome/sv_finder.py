import pysam
import re
import sys
import os
import subprocess
from itertools import groupby
from pavfinder.genome import gapped_align
from pavfinder.genome import split_align
from pavfinder.genome.adjacency import Adjacency
from pavfinder.genome.variant import Variant
from pavfinder.genome.annotate import overlap_pe, parallel_parse_overlaps, annotate_rna_event, annotate_gene_fusion, update_features, get_acen_coords
from pavfinder.genome.alignment import reverse_complement, target_non_canonical
from pavfinder.genome.vcf import VCF


class SVFinder:
    def __init__(self, bam_file, contig_fasta, genome_fasta, out_dir,
                 genome=None, index_dir=None, num_procs=0, skip_simple_repeats=False,
                 cytobands_file=None, acen_buffer=0, debug=False):
        self.bam = pysam.Samfile(bam_file, 'rb')
        self.contig_fasta_file = contig_fasta
        self.contig_fasta = pysam.Fastafile(contig_fasta)
        self.ref_fasta = pysam.Fastafile(genome_fasta)
        self.genome = genome
        self.index_dir = index_dir
        self.num_procs = num_procs
        self.out_dir = out_dir
        self.adjs = []
        self.skip_simple_repeats = skip_simple_repeats
        self.cytobands_file = cytobands_file
        self.acen_buffer = acen_buffer
        self.debug = debug

        self.avg_tlen = None
        self.avg_tlen_normal = None

    def find_adjs(self, min_ctg_cov, max_size=None, min_size=None, ins_as_ins=False,
                  skip_acen=False, check_alt_paths=False, min_ctg_size=0, bad_coords=None, skip_contigs_file=None):
        """Main method to go through the BAM file, extract split and gapped alignments, and calls
        the respective modules to identify adjs"""
        def find_events_in_single_align(align):
            """Implement as sub-function so that small-scale events can be found on split alignments too"""
            adjs = gapped_align.find_adjs(align,
                                          contig_seq,
                                          False,
                                          ins_as_ins=ins_as_ins,
                                          query_fasta=self.contig_fasta,
                                          target_fasta=self.ref_fasta)

            repeats = set()
            for i in range(len(adjs)):
                adj = adjs[i]

                if self.skip_simple_repeats and self.break_region_has_low_complexity(adj.chroms[0], adj.breaks):
                    repeats.add(i)
                    if self.debug:
                        sys.stdout.write("remove contig {} {} potential simple-repeat {}:{}-{}\n".format(adj.contigs[0],
                                                                                                         adj.rearrangement,
                                                                                                         adj.chroms[0],
                                                                                                         adj.breaks[0],
                                                                                                         adj.breaks[1]))
                    continue

            if repeats:
                for i in sorted(repeats, reverse=True):
                    del adjs[i]

            return adjs

        def is_align_in_acen(align, acen):
            """Checks to see if alignment overlaps with acentromeric coordinates
            Args:
                align: alignment (Alignment)
                acen: acentromeric coordinates parsed from UCSC cytobands file (Dictionary) {chrom:(start, end), (start, end)}
            Returns True if overlapped
            """
            s1, e1 = align.tstart, align.tend
            if align.target in acen:
                for (start, end) in acen[align.target]:
                    s2, e2 = int(start) - self.acen_buffer, int(end) + self.acen_buffer
                    if s1 <= e2 and s2 <= e1:
                        return True

            return False

        def create_set(list_file):
            """Creates set from items in a list"""
            subset = set()
            for line in open(list_file, 'r'):
                subset.add(line.strip('\n'))
            return subset

        acen_coords = None
        if skip_acen:
            acen_coords = get_acen_coords(self.cytobands_file)

        skip_contigs = None
        if skip_contigs_file and os.path.exists(skip_contigs_file):
            skip_contigs = create_set(skip_contigs_file)

        all_adjs = []
        for contig, group in groupby(self.bam.fetch(until_eof=True), lambda x: x.qname):
            print('contig', contig)
            alns = list(group)
            contig_seq = self.contig_fasta.fetch(contig)

            if len(contig_seq) < min_ctg_size:
                if self.debug:
                    sys.stdout.write('{}({} bp) less than min contig size {} bp\n'.format(contig,
                                                                                          len(contig_seq),
                                                                                          min_ctg_size))
                continue

            if skip_contigs and contig in skip_contigs:
                if self.debug:
                    sys.stdout.write('{} skipped\n'.format(contig))
                continue

            if len(alns) > 1:
                chimeric_aligns, dubious = split_align.find_chimera(alns,
                                                                    self.bam,
                                                                    min_coverage=min_ctg_cov,
                                                                    check_alt_paths=check_alt_paths,
                                                                    debug=self.debug)
                if chimeric_aligns:
                    if acen_coords:
                        skip = False
                        for align in chimeric_aligns:
                            if acen_coords and is_align_in_acen(align, acen_coords):
                                if self.debug:
                                    sys.stdout.write('skip contig {} because alignment is in centromere {}:{}-{}\n'.format(contig,
                                                                                                                           align.target,
                                                                                                                           align.tstart,
                                                                                                                           align.tend))
                                skip = True
                                break
                        if skip:
                            continue

                    adjs = split_align.find_adjs(
                        chimeric_aligns, contig_seq, dubious=dubious, debug=self.debug)

                    bad = set()
                    for i in range(len(adjs)):
                        adj = adjs[i]
                        # check if homol is simple repeat
                        if adj.homol_seq and adj.homol_seq[0] != '-' and self.is_homol_low_complexity(adj):
                            if self.debug:
                                sys.stdout.write(
                                    "homol_seq is simple-repeat {}:{}\n".format(adj.contigs[0], adj.homol_seq[0]))
                            bad.add(i)

                        # check if event is simple repeat expansions
                        if self.skip_simple_repeats and self.is_novel_sequence_repeat(adj):
                            if self.debug:
                                sys.stdout.write("novel_seq is simple-repeat {}:{}\n".format(adj.contigs[0], adj.novel_seq))
                            bad.add(i)

                        # inversion with size of 1
                        if adj.rearrangement == 'inv' and adj.get_size() <= 1:
                            if self.debug:
                                sys.stdout.write("inversion with unreasonable size {}:{} {}:{}-{}\n".format(adj.contigs[0],
                                                                                                            adj.get_size(),
                                                                                                            adj.chroms[0],
                                                                                                            adj.breaks[0],
                                                                                                            adj.breaks[1]))
                                bad.add(i)

                        if i > 0:
                            if adjs[i].chroms == adjs[i - 1].chroms and\
                               adjs[i].breaks == adjs[i - 1].breaks and\
                               adjs[i].orients == adjs[i].orients and\
                               adjs[i].contig_breaks != adjs[i - 1].contig_breaks:
                                if self.debug:
                                    sys.stdout.write("{} has 2 contig_breaks for same event\n".format(adj.contigs[0]))
                                bad.add(i - 1)
                                bad.add(i)

                    if bad:
                        for i in sorted(bad, reverse=True):
                            del adjs[i]

                    all_adjs.extend(adjs)

                    # capture small-scale events within each chimeric alignment
                    for align in chimeric_aligns:
                        all_adjs.extend(find_events_in_single_align(align))

            best_align = gapped_align.find_single_unique(
                alns, self.bam, debug=self.debug)
            if best_align:
                all_adjs.extend(find_events_in_single_align(best_align))

        merged_adjs = Adjacency.merge(all_adjs)

        # screen out adjacencies that overlap segdups
        if bad_coords is not None and os.path.exists(bad_coords):
            self.screen_by_coordinate(merged_adjs, bad_coords)

        # size filtering
        if max_size is not None or min_size is not None:
            selected = []
            for adj in merged_adjs:
                size = adj.get_size()

                if max_size is not None and min_size is not None:
                    if isinstance(size, int) and size >= min_size and size <= max_size:
                        selected.append(adj)

                elif max_size is not None:
                    if isinstance(size, int) and size <= max_size:
                        selected.append(adj)

                elif min_size is not None:
                    if not isinstance(size, int) or size >= min_size:
                        selected.append(adj)

            return selected
        else:
            return merged_adjs

    def create_variants(self, adjs):
        def track_adjs(used_ids, variants):
            if variants:
                for variant in variants:
                    for adj in variant.adjs:
                        used_ids.add(adj.id)

        """Creates variants from adjacencies"""
        self.variants = []
        adjs_ids_used = set()

        split_events = [adj for adj in adjs if adj.rearrangement not in ('trl', 'ins') and adj.align_types[0] == 'split']
        ins_variants, split_events_remained = Adjacency.extract_interchrom_ins(split_events)
        self.variants.extend(ins_variants)
        track_adjs(adjs_ids_used, ins_variants)

        # special cases for imprecise insertions
        ins_variants, ins_adjs = Adjacency.extract_imprecise_ins([adj for adj in adjs if adj.align_types[0] == 'split' and\
                                                                  adj.rearrangement != 'inv' and\
                                                                  adj.id not in adjs_ids_used
                                                                  ], debug=self.debug)
        self.variants.extend(ins_variants)
        track_adjs(adjs_ids_used, ins_variants)

        # handle inversions
        invs = [adj for adj in adjs if adj.rearrangement == 'inv' and adj.id not in adjs_ids_used]
        inv_variants = Adjacency.group_inversions(invs)
        self.variants.extend(inv_variants)
        track_adjs(adjs_ids_used, inv_variants)

        # convert translocations to insertions
        trls = [adj for adj in adjs if adj.rearrangement == 'trl' and adj.id not in adjs_ids_used]
        ins_variants, trls_remained = Adjacency.extract_interchrom_ins(trls)
        self.variants.extend(ins_variants)
        track_adjs(adjs_ids_used, ins_variants)

        # group reciprocal transcloations
        trls = [adj for adj in adjs if adj.rearrangement == 'trl' and adj.id not in adjs_ids_used]
        reciprocal_trls, trls_remained = Adjacency.group_trls(trls)
        self.variants.extend(reciprocal_trls)
        track_adjs(adjs_ids_used, reciprocal_trls)

        # append remaining non-dubious translocations
        trls = [adj for adj in adjs if adj.rearrangement == 'trl' and adj.id not in adjs_ids_used]
        for trl in trls:
            if not trl.dubious:
                variant = Variant('TRL', [trl])

        for adj in adjs:
            if adj.id not in adjs_ids_used and not adj.dubious:
                self.variants.append(Variant(adj.rearrangement.upper(), [adj]))

    def is_novel_sequence_repeat(self, adj, min_len=3):
        """Check to see if novel sequence is part of a tandem repeat

        Only check if length of novel sequence is least min_len in size

        Args:
            adj: (Adjacency)
            min_len: (int) minimum length of novel sequence before consideration
        """
        contig_seq = self.contig_fasta.fetch(adj.contigs[0]).upper()
        contig_breaks = (adj.contig_breaks[0][0], adj.contig_breaks[0][1] - 2)

        len_seq = 0
        if adj.novel_seq and len(adj.novel_seq) > min_len:
            len_seq = len(adj.novel_seq)
        copies = []

        if len_seq > 0:
            for i in range(len(contig_seq) - len(adj.novel_seq) + 1):
                if i == contig_breaks[0]:
                    continue

                if contig_seq[i:i + len_seq] == adj.novel_seq:
                    # check if it overlaps (or contiguous) with contig_breaks
                    if min(i + len(adj.novel_seq) - 1, contig_breaks[1]) - max(i, contig_breaks[0]) >= -1:
                        copies.append(i)

        return len(copies) > 0

    def break_region_has_low_complexity(self, chrom, breaks, min_units=4, buf=1):
        """Determines if genomic region around deletion/insertion breakpoint has low-complexity sequences

        Window of search = 100bp on either side genomic breakpoint
        True if there is at least "min_units" of repeats that reside in breakpoint region

        Args:
            chrom: (str) chromosome name
            breaks: (List/Tuple) genomic breakpoint (start, end)
            min_units: (int) minimum number of repeat units for returning 'True'
            buf: (int) buffer allowed for considering if repeats reside in breakpoint region

        Returns:
            True if yes, False if no
        """
        r = re.compile(r"(.+?)\1+")
        window = 100
        for i in range(breaks[0] - window, breaks[1] + 1):
            try:
                seq = self.ref_fasta.fetch(chrom, max(0, i), i + 100)
            except BaseException:
                print("can't extract reference sequence for complexity checking {}:{}-{}".format(chrom, i, i + 100))
                continue

            repeats = r.findall(seq)

            if repeats and repeats[0].upper() != 'N':
                m = re.search('^({}){1,}'.format(repeats[0]), seq)
                if m is not None:
                    repeat_start = i + 1
                    repeat_end = i + len(m.group(0))
                    num_units = len(m.group(0)) / len(repeats[0])
                    if num_units >= min_units:
                        if breaks[0] != breaks[1] and breaks[0] + 1 >= repeat_start - \
                                buf and breaks[1] - 1 <= repeat_end + buf:
                            if self.debug:
                                sys.stdout.write('{} {}:{}-{} in low-complexity region {}:{}-{} {}x{}\n'.format('del',
                                                                                                                chrom,
                                                                                                                breaks[0],
                                                                                                                breaks[1],
                                                                                                                chrom,
                                                                                                                repeat_start,
                                                                                                                repeat_end,
                                                                                                                repeats[0].upper(),
                                                                                                                num_units,))

                            return True
                        elif breaks[0] == breaks[1] and breaks[0] >= repeat_start - buf and breaks[0] + 1 <= repeat_end + buf:
                            if self.debug:
                                sys.stdout.write(
                                    '{} {}:{}-{} in low-complexity region {}:{}-{} {}x{}\n'.format('ins',
                                                                                                   chrom,
                                                                                                   breaks[0],
                                                                                                   breaks[1],
                                                                                                   chrom,
                                                                                                   repeat_start,
                                                                                                   repeat_end,
                                                                                                   repeats[0].upper(),
                                                                                                   num_units,))
                            return True
        return False

    def expand_contig_breaks(self, chrom, breaks, contig, contig_breaks, event, debug=False):
        """Expands contig_breaks if repeats reside in breakpoints

        Args:
            chrom: (str) chromosome name
            breaks: (tuple/list) genomic chromosome breaks
            contig: (str) contig ID
            contig_breaks: (tuple/list) contig breaks
            event: (str) 'del' or 'ins'
            debug: (boolean) report debug statements

        Returns:
            tuple of expanded contig breaks
        """
        def extract_repeat(seq):
            """Extracts repeats from given sequence

            Args:
                seq: (str) sequence

            Returns:
                (start, end) of repeat sequence
            """
            repeat = {'start': None, 'end': None}
            if len(seq) == 1:
                repeat['start'] = seq
                repeat['end'] = seq
            else:
                re_start = re.compile(r"^(.+?)\1+")
                re_end = re.compile(r"(.+?)\1+$")

                repeats = re_start.findall(seq)
                if repeats:
                    repeat['start'] = repeats[0]
                repeats = re_end.findall(seq)
                if repeats:
                    repeat['end'] = repeats[0]

            return repeat

        # extract repeats (if any) from 'del' or 'ins'
        contig_breaks_sorted = sorted(contig_breaks)
        pos_strand = True if contig_breaks[0] < contig_breaks[1] else False
        contig_seq = self.contig_fasta.fetch(contig)
        contig_breaks_expanded = [contig_breaks[0], contig_breaks[1]]

        seq = None
        if event == 'del':
            seq = self.ref_fasta.fetch(chrom, breaks[0], breaks[1] - 1)
            if not pos_strand:
                seq = reverse_complement(seq)
        elif event == 'ins':
            seq = contig_seq[contig_breaks_sorted[0] : contig_breaks_sorted[1] - 1]
        if seq is None:
            print(contig, event, 'cannot find seq', seq)
            return None

        repeat = extract_repeat(seq)
        # downstream
        if repeat['end'] is not None:
            seq = repeat['end']
            size = len(seq)
            start = contig_breaks_sorted[1] - 1
            expand = 0
            while start + size <= len(contig_seq):
                next_seq = contig_seq[start: start + size]
                if next_seq.upper() != seq.upper():
                    break
                else:
                    expand += size
                    start += size
            if pos_strand:
                contig_breaks_expanded[1] += expand
            else:
                contig_breaks_expanded[0] += expand
        # upstream
        if repeat['start'] is not None:
            seq = repeat['start']
            size = len(seq)
            start = contig_breaks_sorted[0]
            expand = 0
            while start - size >= 0:
                next_seq = contig_seq[start - size: start]
                if next_seq.upper() != seq.upper():
                    break
                else:
                    expand -= size
                    start -= size
            if pos_strand:
                contig_breaks_expanded[0] += expand
            else:
                contig_breaks_expanded[1] += expand

        if debug and tuple(contig_breaks_expanded) != contig_breaks:
            sys.stdout.write('contig breaks expanded:{} {} -> {}\n'.format(contig, contig_breaks, contig_breaks_expanded))

        return tuple(contig_breaks_expanded)

    def write_subseq(self, adj, out, name_sep):
        """Outputs the sub-sequence of a split alignment to output file

        Args:
            adj: Variant object
            out: Filehandle of output file
            name_sep: (str) Character used to combined various info into query name
        """
        subseqs = adj.extract_subseqs(self.contig_fasta)
        for i in range(len(subseqs)):
            out.write('>{}{}{}{}{}\n{}\n'.format(adj.contigs[0],
                                                 name_sep,
                                                 adj.key(),
                                                 name_sep,
                                                 i,
                                                 subseqs[i]))

    def is_homol_low_complexity(self, adj, min_len=5):
        """Determine if the microhomology sequence of the Adjacency is low-complexity

        Only determines if microhomology sequence is low-complexity if it's at least
        min_len in size

        Args:
            adj: (Adjacency)
            min_len: (int) minimum length of the microhomology sequence in order for
                           complexity to be considered
        """
        if adj.homol_seq[0] is not None and len(adj.homol_seq[0]) >= min_len:
            if self.is_seq_low_complexity(adj.homol_seq[0]):
                return True

        return False

    def is_subseq_low_complexity(self, adj):
        """Determine if any of the 2 sub-sequences is low-complexity

        Args:
            adj: (Adjacency)
        Returns:
            True if any of the 2 sub-sequences is low-complexity
        """
        subseqs = adj.extract_subseqs(self.contig_fasta)

        for subseq in subseqs:
            if self.is_seq_low_complexity(subseq):
                return True

        return False

    def is_seq_low_complexity(self, seq, threshold=0.9):
        """Determines if given sequence is low-complexity

        Low-complexity conditions:
        1. total of any base / length of sequence > threshold
        2. total of bases in same dimer / length of sequence > threshold

        Args:
            seq: (str) sequence
            threshold: (float) minimum fraction of sequence that are same base or dimers
        Returns:
            True if yes, False if no
        """
        bases = ('A', 'G', 'T', 'C')

        for base in bases:
            if float(seq.upper().count(base)) / float(len(seq)) > threshold:
                return True

            # dimer content
            for i in range(len(bases) - 1):
                for j in range(i + 1, len(bases)):
                    if i == j:
                        continue
                    dimer = bases[i] + bases[j]
                    bases_in_dimers = seq.upper().count(dimer) * 2

                    if float(bases_in_dimers) / float(len(seq)) > threshold:
                        return True

        return False

    def screen_realigns(self, use_realigns=False):
        """Realign probe sequences of adjacencies and screen results

        - genome, and index_dir must have been set when object is initialized
        - output is always set to "realign.fa" and "realign.bam"
        - will fail Adjacency if probe sequence can align to single location
        """
        if not self.genome or not self.index_dir:
            return None

        name_sep = '.'
        all_adjs = []
        for variant in self.variants:
            all_adjs.extend(variant.adjs)
        realign_bam_file = Adjacency.realign(all_adjs,
                                             self.out_dir,
                                             probe=True,
                                             contigs_fasta=self.contig_fasta,
                                             name_sep=name_sep,
                                             genome=self.genome,
                                             index_dir=self.index_dir,
                                             num_procs=self.num_procs,
                                             use_realigns=use_realigns,
                                             )
        try:
            bam = pysam.Samfile(realign_bam_file, 'rb')
        except BaseException:
            sys.exit('Error parsing realignment BAM:{}'.format(realign_bam_file))

        # creates mapping from query to variant and Adjacency
        query_to_variant = {}
        for i in range(len(self.variants)):
            for j in range(len(self.variants[i].adjs)):
                adj = self.variants[i].adjs[j]
                query = adj.contigs[0] + name_sep + adj.key()
                query_to_variant[query] = (i, j)

        failed_variants = set()
        for key, group in groupby(bam.fetch(until_eof=True), lambda x: name_sep.join(x.qname.split(name_sep)[:2])):
            alns = list(group)
            variant_idx = query_to_variant[key][0]
            variant = self.variants[variant_idx]
            adj_idx = query_to_variant[key][1]
            adj = variant.adjs[adj_idx]
            adj_aligns = adj.aligns[0]

            indices_to_check = (0, 1)
            if variant.event == 'INS':
                index = None
                for i in (0, 1):
                    if variant.chrom == adj.chroms[i] and (variant.pos[0] == adj.breaks[i] or variant.pos[1] == adj.breaks[i]):
                        index = i
                        break

                if index is not None:
                    indices_to_check = (index,)

            probe_alns = [aln for aln in alns if not aln.qname[-1].isdigit()]
            if not gapped_align.screen_probe_alns(
                    adj_aligns, probe_alns, adj.align_types[0]):
                if self.debug:
                    sys.stdout.write('probe align completely to one location or not aligned with confidence: {}\n'.formatkey)
                failed_variants.add(variant)
                continue

        for failed_var in failed_variants:
            self.variants.remove(failed_var)

    def output(self, reference_url=None, assembly_url=None, insertion_as_breakends=None, header=None):
        """Wrapper function to output Variants and Adjacencies
        Args:
            only_somatic: (boolean) Only outputs somatic variants/adjacencies
            reference_url: (str) reference url to be put in VCF header (optional)
            assembly_url: (str) assembly url to be put in VCF header (optional)
            insertion_as_breakends: (boolean) Output big insertion as breakends
        """
        variants = [variant for variant in self.variants if not variant.filtered_out]

        source = 'NA'
        if header is not None and 'software' in header:
            source = header['software']
        self.output_variants(variants,
                             '{}/variants.vcf'.format(self.out_dir),
                             reference_url=reference_url,
                             assembly_url=assembly_url,
                             insertion_as_breakends=insertion_as_breakends,
                             source=source,
                             )

        adjs = []
        for variant in variants:
            adjs.extend(variant.adjs)

        header = '#{}\n#{} {}'.format(header['software'], header['time'], header['cmd'])
        self.output_adjacencies(adjs,
                                '{}/adjacencies.bedpe'.format(self.out_dir),
                                format='bedpe',
                                header=header)

    def output_variants(self, variants, out_file, reference_url=None, assembly_url=None,
                        insertion_as_breakends=False, source='NA'):
        """Output variants in VCF format
        Args:
            variants: (List) Variants
            out_file: (str) absolute path of output VCF file
            reference_url: (str) reference url to be put in VCF header (optional)
            assembly_url: (str) assembly url to be put in VCF header (optional)
            insertion_as_breakends: (boolean) Output big insertion as breakends
        """
        records = []
        for variant in variants:
            output = variant.as_vcf(self.ref_fasta, insertion_as_sv=not insertion_as_breakends)

            if output is not None and output != '':
                records.extend(output.split('\n'))

        out = open(out_file, 'w')
        out.write('{}\n'.format(VCF.header(source=source, reference_url=reference_url, assembly_url=assembly_url)))
        # sort by chromosome and pos
        records.sort(key=lambda record: (record.split('\t')[0], record.split('\t')[1]))
        for record in records:
            out.write('{}\n'.format(record))
        out.close()

    def output_probes(self, adjs, out_file):
        """Output probe sequences in FASTA format
        Args:
            adjs: (List) Adjacencies
            out_file: (str) absolute path of output FASTA file
        """
        out = open(out_file, 'w')
        for adj in adjs:
            if adj.probes[0] != 'NA':
                out.write('>{} {}\n{}\n'.format(adj.id, len(adj.probes[0]), adj.probes[0]))
        out.close()

    def output_adjacencies(self, adjs, out_file, format, header=None):
        """Output adjacencies in tsv format
        Args:
            adjs: (List) Adjacencies
            out_file: (str) absolute path of output file
            format: (str) either "tab" or "bedpe"
            header: (str) header string
        """
        fn = None
        args = ()
        if format == 'bedpe':
            fn = 'as_bedpe'
        elif format == 'tab':
            fn = 'as_tab'

        if fn is not None:
            out = open(out_file, 'w')
            if header is not None:
                out.write(header + '\n')

            if format == 'tab':
                out.write('{}\n'.format(Adjacency.show_tab_headers()))
            elif format == 'bedpe':
                out.write('{}\n'.format(Adjacency.show_bedpe_headers()))

            for adj in adjs:
                output = getattr(adj, fn)(*args)
                try:
                    out.write('{}\n'.format(output))
                except BaseException:
                    sys.stdout.write("can't output Adjacency")

            out.close()

    def find_support(self, script, bam, min_support, min_overlap,
                     allow_clipped=False, normal_bam=None, min_support_normal=None, min_overlap_normal=None,
                     allow_clipped_normal=False, min_ratio_mapped=None, force=False, debug=False):
        cmd = "python {} {} {} {} --num_procs {}".format(script,
                                                         self.out_dir,
                                                         bam,
                                                         self.contig_fasta_file,
                                                         self.num_procs,)
        if min_support is not None:
            cmd += ' --min_support {}'.format(min_support)
        if allow_clipped:
            cmd += ' --allow_clipped'
            if min_ratio_mapped is not None:
                cmd += ' --support_min_mapped {}'.format(min_ratio_mapped)
        if min_overlap is not None:
            cmd += ' --min_overlap {}'.format(min_overlap)

        # normal
        if normal_bam is not None:
            cmd += ' --normal_bam {}'.format(normal_bam)
        if min_support_normal is not None:
            cmd += ' --min_support_normal {}'.format(min_support_normal)
        if min_overlap_normal is not None:
            cmd += ' --min_overlap_normal {}'.format(min_overlap_normal)
        if allow_clipped_normal:
            cmd += ' --allow_clipped_normal'
        if force:
            cmd += ' --force'
        if debug:
            cmd += ' --debug'

        print(cmd)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()

    def screen_by_coordinate(self, adjs, bad_bed_file):
        def create_bed(outfile):
            out = open(outfile, 'w')
            for adj in adjs:
                out.write('{}\n'.format(adj.as_bed()))
            out.close()

        # creates bed file for all adjacencies
        adjs_bed_file = '{}/adjs.bed'.format(self.out_dir)
        create_bed(adjs_bed_file)
        adjs_bed = BedTool(adjs_bed_file)

        # overlaps adjacencies with bad bed file, keeps breakpoints in Set
        bad_breaks = set()
        regions = BedTool(bad_bed_file)
        overlaps = regions.intersect(adjs_bed)
        for olap in overlaps:
            bad_breaks.add('{}:{}'.format(olap[0], int(olap[1]) + 1))

        # screen out adjacencies
        bad_adj_indices = set()
        for i in range(len(adjs)):
            for j in (0, 1):
                breakpt = '{}:{}'.format(adjs[i].chroms[j], adjs[i].breaks[j])
                if breakpt in bad_breaks:
                    bad_adj_indices.add(i)
                    if self.debug:
                        sys.stdout.write('{} {}:{} {}:{} ({}) overlaps repeat/segdup\n'.format(adjs[i].contigs[0],
                                                                                               adjs[i].chroms[0],
                                                                                               adjs[i].breaks[0],
                                                                                               adjs[i].chroms[1],
                                                                                               adjs[i].breaks[1],
                                                                                               breakpt))
        for i in sorted(bad_adj_indices, reverse=True):
            del adjs[i]

    def filter_variants(self, max_homol=None):
        """Filter out events that are believed to be false positive

        Filtering:
        1. Probe sequence doesn't have 'N'
        2. Novel sequence doesn't have non-AGTC sequence
        3. Will not consider haplotype contigs in assembly
        4. Homology sequence > max_homol

        Args:
            min_spanning: (int) minimum number of spanning reads
            min_flanking: (int) minimum number of flanking pairs when homlogous sequence > max_homol_len
        """

        for variant in self.variants:
            for adj in variant.adjs:
                if 'N' in adj.probes[0].upper():
                    adj.filtered_out = True
                    if self.debug:
                        sys.stdout.write('N_in_probe {} {}\n'.format(adj.contigs, adj.probes[0]))

                if not adj.filtered_out and adj.novel_seq is not None and adj.novel_seq != 'NA' and adj.novel_seq != '-' and\
                        re.search('[^ATGC]', adj.novel_seq, re.IGNORECASE):
                    adj.filtered_out = True
                    if self.debug:
                        sys.stdout.write('non AGTC in novel sequence {} {}\n'.format(adj.contigs, adj.novel_seq))

                if not adj.filtered_out and (target_non_canonical(adj.chroms[0]) or target_non_canonical(adj.chroms[1])):
                    adj.filtered_out = True
                    if self.debug:
                        sys.stdout.write('non canonical chromosome {} {}\n'.format(adj.contigs, adj.chroms))

                if not adj.filtered_out and adj.homol_seq and max_homol is not None:
                    if len(adj.homol_seq[0]) > max_homol:
                        adj.filtered_out = True
                        if self.debug:
                            sys.stdout.write('homolgous sequence length {} too long (>{})\n'.format(len(adj.homol_seq[0]), max_homol))

                if adj.filtered_out:
                    variant.filtered_out = True
