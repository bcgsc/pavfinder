from .alignment import compare_chr, reverse_complement


class Variant:
    def __init__(self, event, adjs, chrom=None, pos=None):
        self.event = event
        self.adjs = adjs

        try:
            self.id = '-'.join([str(adj.id) for adj in adjs])
        except BaseException:
            self.id = None

        self.chrom = chrom
        self.pos = pos
        self.filtered_out = None
        self.filtered_out_normal = None
        self.somatic = False

        self.support = None
        self.support_normal = None
        # final integer counts of read support
        self.final_support = None
        self.final_support_normal = None

    def as_vcf(self, ref_fasta, genomic=True, size_threshold=100, insertion_as_sv=True):
        if len(self.adjs) == 2:
            if self.event == 'INV':
                return self.inversion_as_vcf(ref_fasta)

            elif self.event == 'TRL':
                return self.reciprocal_trl_as_vcf(ref_fasta)

            elif self.event == 'INS' and insertion_as_sv:
                return self.insertion_as_vcf(ref_fasta)

        output = []
        for adj in self.adjs:
            output.append(adj.as_vcf(ref_fasta, size_threshold, genomic))
        return '\n'.join(output)

    def inversion_as_vcf(self, ref_fasta):
        info = {}

        contigs = set()
        for i in range(2):
            for contig in self.adjs[i].contigs:
                contigs.add(contig)
        info['BKPTID'] = ','.join(contigs)

        if self.adjs[0].homol_seq and self.adjs[0].homol_seq[0] != '-':
            info['CIPOS'] = '0,{}'.format(len(self.adjs[0].homol_seq[0]))

        if self.adjs[1].homol_seq and self.adjs[1].homol_seq[0] != '-':
            info['CIEND'] = '0,{}'.format(len(self.adjs[1].homol_seq[0]))

        info['READSUPPORT'] = self.final_support

        if self.somatic:
            info['SOMATIC'] = 'SOMATIC'

        return self.adjs[0].as_sv(ref_fasta, id_ext=self.id, info_ext=info)

    def reciprocal_trl_as_vcf(self, ref_fasta):
        parids = ((self.adjs[1].id + 'a', self.adjs[1].id + 'b'),
                  (self.adjs[0].id + 'a', self.adjs[0].id + 'b'))

        event_id = self.event.upper() + self.id

        info_0 = {'EVENT': (event_id, event_id), 'PARID': parids[0]}
        info_1 = {'EVENT': (event_id, event_id), 'PARID': parids[1]}
        if self.somatic:
            info_0['SOMATIC'] = 'SOMATIC'
            info_1['SOMATIC'] = 'SOMATIC'

        vcf_0 = self.adjs[0].as_breakends(ref_fasta, info_ext=info_0)
        vcf_1 = self.adjs[1].as_breakends(ref_fasta, info_ext=info_1)

        return vcf_0 + '\n' + vcf_1

    def insertion_as_vcf(self, ref_fasta):
        assert self.chrom is not None, 'must specify chrom for big insertion event'
        assert isinstance(self.pos, list) or isinstance(self.pos, tuple), '"pos" for big insertion event must be list or tuple'
        info = {}

        contigs = set()
        for i in range(2):
            for contig in self.adjs[i].contigs:
                contigs.add(contig)
        info['BKPTID'] = ','.join(contigs)
        if self.adjs[0].insertion_size and self.adjs[0].insertion_size > 0:
            info['SVLEN'] = self.adjs[0].insertion_size
        if self.pos[1] != self.pos[0]:
            info['CIPOS'] = '0,%d' % (self.pos[1] - self.pos[0])

        info['READSUPPORT'] = self.final_support

        if self.somatic:
            info['SOMATIC'] = 'SOMATIC'

        return self.adjs[0].as_sv(ref_fasta, chrom_ext=self.chrom, pos_ext=self.pos[0], id_ext=self.id, info_ext=info)

    def check_adjs(self, adjs_passed):
        num_adjs_passed = 0
        for adj in self.adjs:
            if adj.id in adjs_passed:
                num_adjs_passed += 1

        if num_adjs_passed == 2:
            return True
        else:
            return False
