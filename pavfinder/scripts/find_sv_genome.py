import argparse
import sys
import os
import datetime
import pavfinder as pv
from pavfinder.genome.sv_finder import SVFinder


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("c2g", type=str, help="contig-to-genome bam")
    parser.add_argument("contigs_fasta", type=str, help="contigs fasta")
    parser.add_argument("genome_fasta", type=str, help="genome fasta")
    parser.add_argument("outdir", type=str, help="outdir")
    parser.add_argument(
        "--genome",
        dest="genome",
        help="genome prefix of bwa index, use together with index_dir for probe sequence filtering")
    parser.add_argument(
        "--index_dir",
        dest="index_dir",
        help="bwa index directory")
    parser.add_argument("--r2c", type=str, help="reads to genome bam")
    parser.add_argument(
        "--normal_bam",
        type=str,
        help="reads-to-contigs bam file of match normal")
    parser.add_argument(
        "--num_threads",
        help="number of threads/processes. Default:8",
        type=int,
        default=8)
    parser.add_argument("--debug", action="store_true", help="debug mode")
    parser.add_argument(
        "--version", action='version', version='%s %s' %
        (pv.__name__, pv.__version__))

    sv = parser.add_argument_group('variant detection')
    sv.add_argument(
        "--max_homol",
        help="maximum bases of microhomology. Default:25",
        type=int,
        default=25)
    sv.add_argument(
        "--min_ctg_cov",
        help="minimum contig coverage. Default:0.95",
        type=float,
        default=0.95)
    sv.add_argument(
        "--skip_simple_repeats",
        help="skip simple repeats",
        action="store_true",
        default=False)
    sv.add_argument("--max_size", help="maximum size of variant", type=int)
    sv.add_argument("--min_size", help="minimum size of variant", type=int)
    sv.add_argument(
        "--ins_as_ins",
        help="keep small duplications as insertions",
        action="store_true",
        default=False)
    sv.add_argument(
        "--use_realigns",
        help="use existing realignments",
        action="store_true",
        default=False)
    sv.add_argument(
        "--skip_acen",
        help="skip acentromeric regions",
        action="store_true",
        default=False)
    sv.add_argument("--cytobands", help="cytobands file")
    sv.add_argument(
        "--acen_buffer",
        help="buffer added to acen region. Default:100000",
        type=int,
        default=100000)
    sv.add_argument(
        "--check_alt_paths",
        dest="check_alt_paths",
        help="combine primary and secondary alignments to check for alternative path",
        action="store_true",
        default=False)
    sv.add_argument(
        "--min_ctg_size",
        help="minimum contig size. Default:0 (no screening)",
        type=int,
        default=0)
    sv.add_argument(
        "--bad_coords",
        help="BED file for coordinates to screen out e.g. segdups")
    sv.add_argument("--skip_contigs", help="text file of contig names to skip")

    support = parser.add_argument_group('support filtering')
    support.add_argument(
        "--min_support",
        help="minimum read support. Default:2",
        type=int,
        default=2)
    support.add_argument(
        "--min_support_normal",
        help="minimum read support for normal. Default:1",
        type=int,
        default=1)
    support.add_argument(
        "--min_overlap",
        help="minimum breakpoint overlap for identifying read support. Default:4",
        type=int,
        default=4)
    support.add_argument(
        "--min_overlap_normal",
        help="minimum breakpoint overlap for identifying read support. Default:4",
        type=int,
        default=4)
    support.add_argument(
        "--allow_clipped",
        help="allow using clipped reads in gathering read support",
        action="store_true",
        default=False)
    support.add_argument(
        "--allow_clipped_normal",
        help="allow using clipped reads in gathering normal read support",
        action="store_true",
        default=False)
    support.add_argument(
        "--support_min_mapped",
        help="when clipped reads are allowed as read support, minimum ratio of read length mapped Default:0.8",
        type=float,
        default=0.8)
    support.add_argument(
        "--force_support",
        help="overwrite previous read support if any",
        action="store_true",
        default=False)

    output = parser.add_argument_group('output')
    output.add_argument("--reference_url", help="URL of reference")
    output.add_argument("--assembly_url", help="URL of assembly file for VCF")
    output.add_argument(
        "--insertion_as_breakends",
        help="outputs large insertion as breakends in VCF (Default will output as SV)",
        action="store_true",
        default=False)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    sv_finder = SVFinder(args.c2g,
                         args.contigs_fasta,
                         args.genome_fasta,
                         args.outdir,
                         genome=args.genome,
                         index_dir=args.index_dir,
                         num_procs=args.num_threads,
                         skip_simple_repeats=args.skip_simple_repeats,
                         cytobands_file=args.cytobands,
                         acen_buffer=args.acen_buffer,
                         debug=args.debug)

    # discover adjacencies from alignments
    adjs = sv_finder.find_adjs(min_ctg_cov=args.min_ctg_cov,
                               max_size=args.max_size,
                               min_size=args.min_size,
                               ins_as_ins=args.ins_as_ins,
                               skip_acen=args.skip_acen,
                               check_alt_paths=args.check_alt_paths,
                               min_ctg_size=args.min_ctg_size,
                               bad_coords=args.bad_coords,
                               skip_contigs_file=args.skip_contigs)

    # create variants from adjacencies
    sv_finder.create_variants(adjs)
    # filter out variants based on alignments, etc
    sv_finder.filter_variants(max_homol=args.max_homol)

    # realign to increase specificity
    if args.genome and args.index_dir and sv_finder.variants:
        sv_finder.screen_realigns(use_realigns=args.use_realigns)

    # output
    cmd = ' '.join(sys.argv)
    time = datetime.datetime.now().strftime("%Y-%m-%d:%H:%M:%S")
    software = '%s %s' % (pv.__name__, pv.__version__)
    sv_finder.output(reference_url=args.reference_url,
                     assembly_url=args.assembly_url,
                     insertion_as_breakends=args.insertion_as_breakends,
                     header={'cmd': cmd, 'time': time, 'software': software},
                     )

    if args.r2c:
        support_script = '%s/check_support.py' % (os.path.dirname(__file__))
        sv_finder.find_support(
            support_script,
            args.r2c,
            args.min_support,
            args.min_overlap,
            args.allow_clipped,
            args.normal_bam,
            args.min_support_normal,
            args.min_overlap_normal,
            args.allow_clipped_normal,
            args.support_min_mapped,
            force=args.force_support,
            debug=args.debug)


main()
