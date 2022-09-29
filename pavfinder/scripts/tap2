#!/usr/bin/env python
import pavfinder as pv
from ruffus import *
import ruffus.cmdline as cmdline
import subprocess
from collections import defaultdict
import gzip
import os
import sys
import glob
import re
import fileinput
import datetime
import itertools
from configparser import ConfigParser

required_params = {'sv': ['genome_fasta', 'gtf', 'transcripts_fasta', 'genome_index'],
                   'splice': ['genome_fasta', 'gtf']}

def run_cmd(cmd, force=False):
    """Execute command"""
    process = subprocess.Popen(cmd,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE,
                               shell = True)

    stdout, stderr = process.communicate()

    if process.returncode != 0 and not force:
        raise Exception("Failed to run '{}'\n{}{}Non-zero exit status {}".format(cmd,
                                                                                 stdout.decode('utf-8'),
                                                                                 stderr.decode('utf-8'),
                                                                                 process.returncode))

    return stdout.decode('utf-8'), stderr.decode('utf-8')
    
def format_read_pairs(fqs=None, list_file=None):
    """Formulates fastq pairs from either fastq parameter or file that lists fastq files"""
    fastqs = []
    if list_file is not None:
        with open(list_file, 'r') as ff:
            for line in ff:
                fastqs.append(line.rstrip('\n'))
                
    elif fastqs is not None:
        fastqs = fqs
            
    fastqs1 = sorted([fq for fq in fastqs if '1.' in fq])
    fastqs2 = sorted([fq for fq in fastqs if '2.' in fq])
    
    if len(fastqs1) != len(fastqs2):
        raise Exception('input fastqs not paired')
    
    return fastqs1, fastqs2
                
def get_version(exe):
    """Gets software version for formulating directory names"""
    version = None
    cmd = None
    if exe == 'bbt':
        cmd = 'biobloommicategorizer --version'
        re.compile(r'(\d+\.\d+\.\d+[a-z])')
    elif exe == 'rnabloom':
        cmd = 'rnabloom -v'
    elif exe == 'pv':
        version = pv.__version__

    if cmd is not None:
        stdout, stderr = run_cmd(cmd)
        # bbt output to stderr, rnabloom to stdout
        match = None
        if stderr:
            match = re.search(r'(\d+\.\d+\.\d+[a-z]?)', stderr)
        if not match and stdout:
            match = re.search(r'(\d+\.\d+\.\d+[a-z]?)', stdout)
        if match:
            version = match.group(1)

    return version

def get_params(args, params_file):
    """Extracts 'alignments' and 'annotations' parameters from either args or params_file"""
    def override(section, names):
        for name in names:
            value = getattr(args, name)
            if not section in params:
                params[section] = {}
            if not name in params[section] or value is not None:
                params[section][name] = value

    params = {}
    if params_file is not None and os.path.exists(params_file):
        cfg = ConfigParser()
        cfg.read(params_file)
        for section in cfg.sections():
            params[section] = {}
            for option in cfg.options(section):
                value = cfg.get(section, option)
                if value.lower() in ('true', 'false'):
                    if value.lower() == 'true':
                        value = True
                    else:
                        value = False
                elif ' ' in value:
                    value = value.split()
                elif value.isspace():
                    value = None
                params[section][option] = value

    override('alignments', ('genome_index', 'transcripts_fasta', 'sort_mem'))
    override('annotations', ('genome_fasta', 'gtf', 'suppl_annot'))

    return params

def check_params(args, params):
    """Check if required 'alignments' and 'annotations' parameters are specified"""
    error = None
    if not args.only_assembly:
        if not 'alignments' in params:
            error = 'no alignments parameters provided'
        elif not 'genome_index' in params['alignments'] or params['alignments']['genome_index'] is None or\
             type(params['alignments']['genome_index']) is not list or len(params['alignments']['genome_index']) != 2:
            error = 'no proper GMAP index provided'
        elif not 'genome_fasta' in params['annotations'] or params['annotations']['genome_fasta'] is None:
            error = 'no genome fasta provided'
        elif not 'gtf' in params['annotations'] or params['annotations']['gtf'] is None:
            error = 'no gtf provided'

        if not args.only_splicing:
            if not 'transcripts_fasta' in params['alignments'] or params['alignments']['transcripts_fasta'] is None:
                error = 'no transcripts fasta provided'

    if error is not None:
        sys.exit('TAP ABORTED: {}'.format(error.upper()))

parser = cmdline.get_argparse(description='TAP pipeline')
parser.add_argument('sample', type=str, help='sample name')
parser.add_argument('outdir', type=str, help='output directory')
parser.add_argument('--bf', type=str, help='path to bloomfilter')
parser.add_argument('--fq', type=str, nargs='+', help='input gzipped fastqs')
parser.add_argument('--fq_list', type=str, help='text file of input fastq paths')
parser.add_argument('--nprocs', type=int, default=32, help='number of threads/processes. Default=32')
parser.add_argument('--remove_fq', action='store_true', help='remove intermediate fastqs')
parser.add_argument('--only_assembly', action='store_true')
parser.add_argument('--only_sv', action='store_true')
parser.add_argument('--only_splicing', action='store_true')
parser.add_argument('--genome_bam', type=str, help='genome bam(for detecting splice-site variants)')
parser.add_argument('--params', type=str, help='parameters file')
assembly = parser.add_argument_group('assembly')
assembly.add_argument('--readlen', type=int, help='read length')
alignments = parser.add_argument_group('alignments')
alignments.add_argument('--genome_index', type=str, nargs=2, help='gmap index')
alignments.add_argument('--transcripts_fasta', type=str, help='bwa index of transcript sequences')
alignments.add_argument('--sort_mem', type=str, help='samtools sort memory. Default:5G', default='5G')
annotations = parser.add_argument_group('annotations')
annotations.add_argument('--genome_fasta', type=str, help='genome fasta')
annotations.add_argument('--gtf', type=str, help='gtf')
annotations.add_argument("--suppl_annot", type=str, nargs="+", help="supplementary annotation file(s) for checking novel splice events")

args = parser.parse_args()
params = get_params(args, args.params)
check_params(args, params)

logs_dir = args.outdir + '/logs'
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

log_file = '{}/log.{}.txt'.format(logs_dir, datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))
logger, logging_mutex = cmdline.setup_logging(__name__,
                                              log_file,
                                              args.verbose)
print('log_file:', log_file)

cmdline.run(args)

read_pairs = []
if args.fq:
    read_pairs = format_read_pairs(fqs=args.fq)
elif args.fq_list:
    read_pairs = format_read_pairs(list_file=args.fq_list)

history_file = '{}/.ruffus_history.sqlite'.format(args.outdir)
bbt_outdir = '{}/bbt_{}'.format(args.outdir, get_version('bbt'))
assembly_outdir = '{}/rnabloom_{}'.format(args.outdir, get_version('rnabloom'))
pv_outdir = '{}/pv_{}'.format(args.outdir, get_version('pv'))
bbt_prefix = bbt_outdir + '/' + args.sample

# for determining how many procs/threads to give to each analysis
num_analysis = 2
if args.only_sv or args.only_splicing:
    num_analysis = 1
elif args.only_assembly:
    num_analysis = 0

assembly_input = []

@follows(mkdir(bbt_outdir))
@transform([read_pairs],
           formatter('_1.fastq.gz$'),
           bbt_outdir + '/' + args.sample + '_reads.fq',
           bbt_prefix,
           args.bf,
           args.nprocs)
def classify(paired_fqs, bbt_output, prefix, bf, nthreads):
    """Classifies reads by gene"""
    if bf:
        if len(paired_fqs[0]) == 1:
            input_fq1 = paired_fqs[0][0]
        else:
            input_fq1 = '<(zcat {})'.format(' '.join(paired_fqs[0]))
        if len(paired_fqs[1]) == 1:
            input_fq2 = paired_fqs[1][0]
        else:
            input_fq2 = '<(zcat {})'.format(' '.join(paired_fqs[1]))

        cmd = 'biobloommicategorizer -t {} -f {} --fq -p {} -e {} {}'.format(nthreads,
                                                                             bf,
                                                                             prefix,
                                                                             input_fq1,
                                                                             input_fq2)
        stdout, stderr = run_cmd('/bin/bash -c "{}"'.format(cmd))

        if stdout:
            with open('{}/{}_reads.fq'.format(bbt_outdir, args.sample), 'w') as out:
                out.write('{}'.format(stdout))

    else:
        cmd = 'touch {}'.format(bbt_output)
        run_cmd('/bin/bash -c "{}"'.format(cmd))
    
def format_read_pairs_for_abyss(lines):
    """format reads to '/1' and '/2' for abyss"""
    read_name1 = lines[0].split()[0]
    read_name2 = lines[4].split()[0]

    if not (read_name1[-2:] == '/1' and read_name2[-2:] == '/2'):
        if read_name1 == read_name2:
            lines[0] = lines[0].replace(read_name1, read_name1 + '/1', 1)
            lines[4] = lines[4].replace(read_name2, read_name2 + '/2', 1)

        elif read_name1[-1] == '1' and read_name2[-1] == '2' and\
             not read_name1[-2].isalnum() and not read_name2[-2].isalnum():
            lines[0] = lines[0].replace(read_name1, read_name1[:-2] + '/' + read_name1[-1:], 1)
            lines[4] = lines[4].replace(read_name2, read_name2[:-2] + '/' + read_name2[-1:], 1)

    return lines

@split(classify,
       bbt_outdir + '/*.fastq*')
def split_input(bbt_fastq, split_fastqs, genes=None):
    """Divides single BBT fastq into single files for each gene""" 
    seqs1 = defaultdict(list)
    seqs2 = defaultdict(list)

    if args.bf:
        with open(bbt_fastq[0], 'r') as fq:
            for lines in itertools.zip_longest(*[fq]*8):
                header_cols = lines[0].rstrip().split('\t')
                if len(header_cols) == 3:
                    targets = [hit.split(',')[0].replace('.fa', '') for hit in header_cols[2].split(';')]
                    lines = format_read_pairs_for_abyss(list(lines))
                    for target in targets:
                        seqs1[target].extend(lines[:4])
                        seqs2[target].extend(lines[-4:])

        split_fastqs = output_split_pairs(seqs1, seqs2, bbt_outdir)

    else:
        for i in range(len(read_pairs[0])):
            source = os.path.abspath(read_pairs[0][i])
            if len(read_pairs[0]) == 1:
                target = '{}/{}_{}.fastq'.format(bbt_outdir, args.sample, 1)
            else:
                target = '{}/{}_{}-{}.fastq'.format(bbt_outdir, args.sample, 1, i + 1)
            if os.path.splitext(source)[1] == '.gz':
                target += '.gz'

            if not os.path.exists(target):
                os.symlink(source, target)

            source = os.path.abspath(read_pairs[1][i])
            if len(read_pairs[1]) == 1:
                target = '{}/{}_{}.fastq'.format(bbt_outdir, args.sample, 2)
            else:
                target = '{}/{}_{}-{}.fastq'.format(bbt_outdir, args.sample, 2, i + 1)
            if os.path.splitext(source)[1] == '.gz':
                target += '.gz'

            if not os.path.exists(target):
                os.symlink(source, target)

def output_split_pairs(seqs1, seqs2, outdir):
    """Outputs read sequences into FASTA files"""
    fqs = []
    for target in sorted(seqs1.keys()):
        fq1 = '{}/{}_1.fastq'.format(outdir, target)
        fq2 = '{}/{}_2.fastq'.format(outdir, target)
        with open(fq1, 'w') as out1:
            out1.writelines(seqs1[target])
                
        with open(fq2, 'w') as out2:
            out2.writelines(seqs2[target])

        fqs.append(fq1)
        fqs.append(fq2)
        
    return fqs

@collate(split_input,
         formatter(".+/(.+?)_[12]-*\d*.fastq.*"),
         assembly_outdir + "/{1[0]}/{1[0]}.transcripts.fa",
         args.readlen)
def assemble_single_gene(fastqs, gene_fasta, readlen):
    outdir = os.path.dirname(gene_fasta)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    prefix = os.path.basename(outdir)
    
    fastqs_sorted = sorted(fastqs)
    
    cmd = 'rnabloom -length {} -l {} -r {} -o {} -n {} -prefix {} --revcomp-right -tiplength 5 -extend -chimera'.format(readlen,
                                                                                          fastqs_sorted[0],
                                                                                          fastqs_sorted[1],
                                                                                          outdir,
                                                                                          prefix,
                                                                                          prefix + '.',
                                                                                          )

    if 'assembly' in params:
        aps = params['assembly']
        cmd += '{} '.format(' '.join(['-{} {}'.format(name, aps[name]) for name in aps.keys() if type(aps[name]) is not bool]))
        cmd += '{} '.format(' '.join(['-{}'.format(name) for name in aps.keys() if type(aps[name]) is bool and aps[name] == True]))

    # rna-bloom will crash if fastq is emtpy
    if os.stat(fastqs_sorted[0]).st_size == 0:
        flag = '{}/TRANSCRIPTS.DONE'.format(outdir)
        cmd = 'touch {} && touch {}'.format(gene_fasta, flag)
        run_cmd(cmd, force=True)
    else:
        stdout_str, stderr_str = run_cmd(cmd, force=True)
        
def bbt_cleanup():
    """Removes all fastq files generated by BBT step r2c is done"""
    if args.remove_fq:
        temp_files = []
        for ff in glob.glob('{}/*.fastq'.format(bbt_outdir)):
            temp_files.append(ff)
        for ff in glob.glob('{}/*.fq'.format(bbt_outdir)):
            temp_files.append(ff)

        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)

@active_if(not args.only_assembly)
@transform(assemble_single_gene,
           formatter(".+transcripts.fa"),
           "{0[0]}.bwt")
def r2c_bwa_index(fasta, index):
    """BWA index each gene FASTA for r2c"""
    if os.path.getsize(fasta) > 0:
        cmd = 'bwa index {}'.format(fasta)

    else:
        cmd = 'touch {}'.format(index)
        
    run_cmd(cmd)
    
@active_if(not args.only_assembly)
@transform(r2c_bwa_index,
           formatter(),
           "{path[0]}/r2c.bam")
def r2c(index, r2c_bam):
    gene = list(filter(None, index.split(os.sep)))[-2]
    reads1 = '{}/{}_1.fastq'.format(bbt_outdir, gene)
    reads2 = '{}/{}_2.fastq'.format(bbt_outdir, gene)
    r2c_bam = os.path.dirname(index) + '/r2c.bam'
    
    cmd = 'bwa mem {} {} {} | samtools view -bhS - -o {}'.format(os.path.splitext(index)[0],
                                                                 reads1,
                                                                 reads2,
                                                                 r2c_bam)

    run_cmd('/bin/bash -c "{}"'.format(cmd))

@merge(assemble_single_gene,
       '{}/{}.fa'.format(assembly_outdir, args.sample))
def concat_fasta(gene_fastas, single_merged_fasta):
    """Concatenates every gene assembly into single fasta file"""
    if args.bf:
        fin = fileinput.input(gene_fastas)
        with open(single_merged_fasta, 'w') as fout:
            for line in fin:
                fout.write(line)
        fin.close()
    else:
        source = os.path.relpath(gene_fastas[0], os.path.dirname(single_merged_fasta))
        os.symlink(source, single_merged_fasta)

    cmd = 'samtools faidx {}'.format(single_merged_fasta)
    run_cmd('/bin/bash -c "{}"'.format(cmd))

@active_if(not args.only_assembly)
@merge(r2c,
       '{}/r2c.bam'.format(assembly_outdir),
       )
def r2c_concat(r2c_bams, r2c_cat_bam):
    """Concatenates all r2c bam files into single file"""
    header_file = '{}/r2c_cat.header'.format(assembly_outdir)
    sam_file = '{}/r2c_cat.sam'.format(assembly_outdir)
    if len(r2c_bams) > 1:
        with open(header_file, 'w') as header, open(sam_file, 'w') as sam:
            for i in range(len(r2c_bams)):
                if os.path.getsize(r2c_bams[i]) == 0:
                        continue
                cmd = subprocess.Popen('samtools view -H {}'.format(r2c_bams[i]), shell=True, stdout=subprocess.PIPE)
                for line in cmd.stdout:
                    line = line.decode('utf-8')
                    if i == 0 and line[:3] == '@HD':
                        header.write(line)
                    elif line[:3] == '@SQ':
                        header.write(line)
    
                cmd = subprocess.Popen('samtools view {}'.format(r2c_bams[i]), shell=True, stdout=subprocess.PIPE)
                for line in cmd.stdout:
                    line = line.decode('utf-8')
                    sam.write(line)

        cmd = 'cat {} {} | samtools view -Su - | samtools sort -m {} - -o {}'.format(header_file,
                                                                                     sam_file,
                                                                                     params['alignments']['sort_mem'],
                                                                                     r2c_cat_bam)

        run_cmd('/bin/bash -c "{}"'.format(cmd))

    elif len(r2c_bams) == 1:
        source = os.path.relpath(r2c_bams[0], os.path.dirname(r2c_cat_bam))
        cmd = 'samtools sort -m {} {} -o {}'.format(params['alignments']['sort_mem'],
                                                    r2c_bams[0],
                                                    r2c_cat_bam)
    
        run_cmd('/bin/bash -c "{}"'.format(cmd))        

def r2c_cleanup():
    """Removes individual r2c indices and intermediate files"""
    temp_files = ['{}/r2c_cat.header'.format(assembly_outdir), '{}/r2c_cat.sam'.format(assembly_outdir)]
    for ff in glob.glob('{}/*/*-merged.fa.*'.format(assembly_outdir)):
        temp_files.append(ff)

    if args.bf:
        r2c_bams = []
        for ff in glob.glob('{}/*/r2c.bam'.format(assembly_outdir)):
            r2c_bams.append(ff)
        if len(r2c_bams) > 1:
            temp_files.extend(r2c_bams)

    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)

@active_if(not args.only_assembly)
@posttask(r2c_cleanup)
@posttask(bbt_cleanup)
@transform(r2c_concat,
           suffix('.bam'),
           '.bam.bai')
def r2c_index_concat(r2c_cat_sorted_bam, r2c_cat_sorted_bam_index):
    """Samtools index sorted r2c bam"""
    cmd = 'samtools index {}'.format(r2c_cat_sorted_bam)
    run_cmd(cmd)
           
@active_if(not args.only_assembly)
@transform(concat_fasta,
           formatter(".+.fa$"),
           "{path[0]}/c2g.bam",
           args.genome_index,
           args.nprocs)
def c2g(contigs_fasta, c2g_bam, genome_index, nthreads):
    """Aligns contigs to genome using GMAP"""
    if os.path.getsize(contigs_fasta) > 0:
        cmd = 'gmap -D {} -d {} {} -t {} -f samse -n 0 -x 10 | samtools view -bhS - -o {}'.format(params['alignments']['genome_index'][0],
                                                                                                  params['alignments']['genome_index'][1],
                                                                                                  contigs_fasta,
                                                                                                  nthreads,
                                                                                                  c2g_bam)
    else:
        cmd = 'touch {}'.format(c2g_bam)

    run_cmd('/bin/bash -c "{}"'.format(cmd))
    
@active_if(not args.only_splicing and not args.only_assembly)
@transform(concat_fasta,
           formatter(".+.fa$"),
           "{path[0]}/c2t.bam",
           args.nprocs)
def c2t(contigs_fasta, c2t_bam, nthreads):
    """Aligns contigs to transcripts using BWA mem"""
    if os.path.getsize(contigs_fasta) > 0:
        cmd = 'bwa mem -t {} {} {} | samtools view -bhS - -o {}'.format(nthreads,
                                                                        params['alignments']['transcripts_fasta'],
                                                                        contigs_fasta,
                                                                        c2t_bam)
    else:
        cmd = 'touch {}'.format(c2t_bam)

    run_cmd('/bin/bash -c "{}"'.format(cmd))
    
@active_if(not args.only_splicing and not args.only_assembly)
@follows(mkdir(pv_outdir))
@merge([concat_fasta, c2g, c2t, r2c_index_concat],
       pv_outdir + '/sv.bedpe',
       args.nprocs,
       )
def find_sv(inputs, events_output, nprocs):
    """Finds structural variants using PAVFinder_transcriptome"""
    merged_fasta, c2g_bam, c2t_bam, r2c_index = inputs

    if num_analysis > 0:
        nprocs /= num_analysis

    if os.path.getsize(merged_fasta) > 0:
        cmd = 'find_sv_transcriptome.py --gbam {} --tbam {} --transcripts_fasta {} --genome_index {} --r2c {} --nproc {}'.format(c2g_bam,
                                                                                                                                 c2t_bam,
                                                                                                                                 params['alignments']['transcripts_fasta'],
                                                                                                                                 ' '.join(params['alignments']['genome_index']),
                                                                                                                                 os.path.splitext(r2c_index)[0],
                                                                                                                                 int(nprocs))

        if 'sv' in params:
            sparams = params['sv']
            cmd += ' {}'.format(' '.join(['--{} {}'.format(name, sparams[name]) for name in sparams.keys() if type(sparams[name]) is not bool]))
            cmd += ' {}'.format(' '.join(['--{}'.format(name) for name in sparams.keys() if type(sparams[name]) is bool and sparams[name] == True]))

        cmd += ' {}'.format(' '.join([merged_fasta, params['annotations']['gtf'], params['annotations']['genome_fasta'], os.path.dirname(events_output)]))

    else:
        cmd = 'touch {}'.format(events_output)

    run_cmd(cmd)

@active_if(not args.only_sv and not args.only_assembly)
@follows(mkdir(pv_outdir))
@merge([concat_fasta, c2g, r2c_index_concat],
       ['{}/{}'.format(pv_outdir, ff) for ff in ('novel_splicing.bedpe', 'mappings.tsv', 'junctions.bed')],
       args.gtf,
       args.genome_fasta,
       args.nprocs,
       args.genome_bam,
       )
def map_splicing(inputs, outputs, gtf, genome_fasta, nprocs, genome_bam):
    """Finds splice_variants, generates coverage and junctions files using PAVFinder_transcriptome"""
    merged_fasta, c2g_bam, r2c_index = inputs

    if num_analysis > 0:
        nprocs /= num_analysis

    if os.path.getsize(merged_fasta) > 0:
        cmd = 'map_splice.py {} {} {} {} {} --r2c {} --nproc {}'.format(c2g_bam,
                                                                        merged_fasta,
                                                                        params['annotations']['gtf'],
                                                                        params['annotations']['genome_fasta'],
                                                                        pv_outdir,
                                                                        os.path.splitext(r2c_index)[0],
                                                                        int(nprocs),
                                                                        )

        if 'splicing' in params:
            sparams = params['splicing']
            cmd += ' {}'.format(' '.join(['--{} {}'.format(name, sparams[name]) for name in sparams.keys() if type(sparams[name]) is not bool]))
            cmd += ' {}'.format(' '.join(['--{}'.format(name) for name in sparams.keys() if type(sparams[name]) is bool and sparams[name] == True]))

        if 'suppl_annot' in params['annotations'] and params['annotations']['suppl_annot'] and\
           params['annotations']['suppl_annot'] is not None:
            if type(params['annotations']['suppl_annot']) is str and not params['annotations']['suppl_annot'].isspace():
                cmd += ' --suppl_annot {}'.format(params['annotations']['suppl_annot'])
            elif type(params['annotations']['suppl_annot']) is list or type(params['annotations']['suppl_annot']) is tuple:
                cmd += ' --suppl_annot {}'.format(' '.join(params['annotations']['suppl_annot']))

        if genome_bam:
            cmd += ' --genome_bam {}'.format(genome_bam)

        run_cmd(cmd)
    else:
        for output in outputs:
            run_cmd('touch {}'.format(output))

pipeline_printout(sys.stdout, verbose=3)
pipeline_run(verbose=3, multiprocess=args.nprocs, logger=logger, history_file=history_file)
