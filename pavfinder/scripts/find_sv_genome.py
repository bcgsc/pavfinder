from optparse import OptionParser
import sys
import os
import datetime
import pavfinder as pv
from pavfinder.genome.sv_finder import SVFinder
					    
def main(args, options):
    sv_finder = SVFinder(*args, 
                         genome=options.genome, 
                         index_dir=options.index_dir, 
                         num_procs=options.num_threads, 
                         skip_simple_repeats=options.skip_simple_repeats,
                         cytobands_file=options.cytobands_file,
                         acen_buffer=options.acen_buffer,
                         debug=options.debug)
        
    # discover adjacencies from alignments
    adjs = sv_finder.find_adjs(min_ctg_cov=options.min_ctg_cov, 
                               max_size=options.max_size, 
                               min_size=options.min_size, 
                               ins_as_ins=options.ins_as_ins,
                               skip_acen=options.skip_acen,
                               check_alt_paths=options.check_alt_paths,
                               min_ctg_size=options.min_ctg_size,
                               bad_coords=options.bad_coords,
                               skip_contigs_file=options.skip_contigs)
    
    # create variants from adjacencies
    sv_finder.create_variants(adjs)
    # filter out variants based on alignments, etc
    sv_finder.filter_variants(max_homol=options.max_homol)
        
    # realign to increase specificity
    if options.genome and options.index_dir and sv_finder.variants:
	sv_finder.screen_realigns(use_realigns=options.use_realigns)
            
    # output
    cmd = ' '.join(sys.argv)
    time = datetime.datetime.now().strftime("%Y-%m-%d:%H:%M:%S")
    software = '%s %s' % (pv.__name__, pv.__version__)
    sv_finder.output(reference_url=options.reference_url,
                     assembly_url=options.assembly_url,
                     insertion_as_breakends=options.insertion_as_breakends,
                     header={'cmd': cmd, 'time':time, 'software':software},
                     )
    
    if options.r2c_bam_file:
	support_script = '%s/check_support.py' % (os.path.dirname(__file__))
	sv_finder.find_support(support_script, options.r2c_bam_file,
	                       options.min_support, options.min_overlap, options.allow_clipped,
	                       options.normal_bam, options.min_support_normal, options.min_overlap_normal,
	                       options.allow_clipped_normal, options.support_min_mapped, force=options.force_support, debug=options.debug)

if __name__ == '__main__':
    usage = "Usage: %prog c2g_bam contig_fasta(indexed) genome_file(indexed) out_dir"
    parser = OptionParser(usage=usage, version=pv.__version__)
    
    parser.add_option("-b", "--r2c_bam", dest="r2c_bam_file", help="reads-to-contigs bam file")
    parser.add_option("-g", "--genome", dest="genome", help="genome")
    parser.add_option("-G", "--index_dir", dest="index_dir", help="genome index directory")
    parser.add_option("-t", "--num_threads", dest="num_threads", help="number of threads. Default:8", type='int', default=8)
    parser.add_option("--debug", dest="debug", help="debug mode", action="store_true", default=False)
    parser.add_option("--max_homol", dest="max_homol", help="maximum bases of microhomology. Default:25", type = "int", default=25)
    parser.add_option("--min_ctg_cov", dest="min_ctg_cov", help="minimum contig coverage. Default:0.95", type='float', default=0.95)
    parser.add_option("--min_support", dest="min_support", help="minimum read support. Default:2", type='int', default=2)
    parser.add_option("--min_support_normal", dest="min_support_normal", help="minimum read support for normal. Default:1", type='int', default=1)
    parser.add_option("--min_overlap", dest="min_overlap", help="minimum breakpoint overlap for identifying read support. Default:4", type='int', default=4)
    parser.add_option("--min_overlap_normal", dest="min_overlap_normal", help="minimum breakpoint overlap for identifying read support. Default:4", type='int', default=4)
    parser.add_option("--allow_clipped", dest="allow_clipped", help="allow using clipped reads in gathering read support", action="store_true", default=False)
    parser.add_option("--allow_clipped_normal", dest="allow_clipped_normal", help="allow using clipped reads in gathering normal read support", action="store_true", default=False)
    parser.add_option("--support_min_mapped", dest="support_min_mapped", help="when clipped reads are allowed as read support, minimum ratio of read length mapped Default:0.8", type='float', default=0.8)
    parser.add_option("--normal_bam", dest="normal_bam", help="reads-to-contigs bam file of match normal")
    parser.add_option("--reference_url", dest="reference_url", help="URL of reference")
    parser.add_option("--assembly_url", dest="assembly_url", help="URL of assembly file for VCF")
    parser.add_option("--skip_simple_repeats", dest="skip_simple_repeats", help="skip simple repeats", action="store_true", default=False)
    parser.add_option("--insertion_as_breakends", dest="insertion_as_breakends", help="outputs large insertion as breakends in VCF (Default will output as SV)", action="store_true", default=False)
    parser.add_option("--max_size", dest="max_size", help="maximum size of variant", type='int')
    parser.add_option("--min_size", dest="min_size", help="minimum size of variant", type='int')
    parser.add_option("--ins_as_ins", dest="ins_as_ins", help="keep small duplications as insertions", action="store_true", default=False)
    parser.add_option("--use_realigns", dest="use_realigns", help="use existing realignments", action="store_true", default=False)
    parser.add_option("--skip_acen", dest="skip_acen", help="skip acentromeric regions", action="store_true", default=False)
    parser.add_option("--cytobands", dest="cytobands_file", help="cytobands file")
    parser.add_option("--acen_buffer", dest="acen_buffer", help="buffer added to acen region. Default:100000", type="int", default=100000)
    parser.add_option("--check_alt_paths", dest="check_alt_paths", help="combine primary and secondary alignments to check for alternative path", action="store_true", default=False)
    parser.add_option("--min_ctg_size", dest="min_ctg_size", help="minimum contig size. Default:0 (no screening)", type="int", default=0)
    parser.add_option("--bad_coords", dest="bad_coords", help="BED file for coordinates to screen out e.g. segdups")
    parser.add_option("--skip_contigs", dest="skip_contigs", help="text file of contig names to skip")
    parser.add_option("--force_support", dest="force_support", help="overwrite previous read support if any", action="store_true", default=False)
    
    (options, args) = parser.parse_args()
    if len(args) == 4:
        main(args, options)


