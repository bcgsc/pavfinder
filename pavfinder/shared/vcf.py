from collections import OrderedDict
import time

class VCF:
    # 8 mandatory fields. 'FORMAT' is optional
    fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
    info_fields = ('SVTYPE', 'EVENTTYPE', 'PARID', 'MATEID', 'END', 'SVLEN', 'CIPOS', 'HOMLEN', 'HOMSEQ', 'BKPTID', 'EVENT', 'SPANNING_READS', 'FLANKING_PAIRS', 'SOMATIC',
                   'REPEAT_SEQ', 'REPEAT_NUM', 'REPEAT_NUM_CHANGE')
    
    version = '4.2'
    
    # INFO key=ID, value=(Number, Type, Description)
    meta_info = OrderedDict()
    meta_info['SVTYPE'] = (1, 'String', 'Type of structural variant')
    meta_info['EVENTTYPE'] = (1, 'String', 'Type of structural variant when breakend notation is used')
    meta_info['PARID'] = (1, 'String', 'ID of partner breakend')
    meta_info['MATEID'] = ('.', 'String', 'ID of mate breakends')
    meta_info['END'] = (1, 'Integer', 'End position of the variant described in this record')
    meta_info['SVLEN'] = ('.', 'Integer', 'Difference in length between REF and ALT alleles')
    meta_info['CIPOS'] = (2, 'Integer', 'Confidence interval around POS for imprecise variants')
    meta_info['HOMLEN'] = ('.', 'Integer', 'Length of base pair identical micro-homology at event breakpoints')
    meta_info['HOMSEQ'] = ('.', 'String', 'Sequence of base pair identical micro-homology at event breakpoints')
    meta_info['BKPTID'] = ('.', 'String', 'ID of the assembled alternate allele in the assembly file')
    meta_info['EVENT'] = ('1', 'String', 'ID of event associated to breakend')
    meta_info['SPANNING_READS'] = (1, 'Integer', 'Number of reads spanning breakpoint')
    meta_info['FLANKING_PAIRS'] = (1, 'Integer', 'Number of unique read pairs flanking breakpoint')
    meta_info['SOMATIC'] = (0, 'Flag', 'Somatic')
    meta_info['REPEAT_SEQ'] = (1, 'Integer', 'Repeat sequence in tandem duplication')
    meta_info['REPEAT_NUM'] = (1, 'Integer', 'Number of novel repeats in tandem duplication')
    meta_info['REPEAT_NUM_CHANGE'] = (1, 'String', 'Change of repeat number in tandem duplication')
    
    # ALT  key=ID, value=Description
    meta_alt = OrderedDict()
    meta_alt['DEL'] = 'Deletion'
    meta_alt['INS'] = 'Insertion of novel sequence'
    meta_alt['INV'] = 'Inversion'
    meta_alt['DUP:TANDEM'] = 'Tandem Duplication'
    
    def __init__(self, chrom, pos, id, ref, alt, info=None):
        self.chrom, self.pos, self.id, self.ref, self.alt = chrom, pos, id, ref, alt
	        
        if info is None:
            self.info = {}
            for key in self.meta_info.keys():
                self.info[key] = None
        else:
            self.info = info
            
    def output(self): 
        """Outputs an invididual VCF record
        Returns:
            string of VCF record
        """  
        data = []
        for field in self.fields:
            if not hasattr(self, field.lower()):
                data.append('.')
            
            elif field == 'INFO':
                info = []
                for item in self.meta_info.keys():
                    if self.info.has_key(item):
                        if self.meta_info[item][1] == 'Flag':
                            info.append(item)
                        else:
                            info.append('%s=%s' % (item, self.info[item]))
                if info:
                    data.append(';'.join(info))
                else:
                    data.append('.')
                        
            else:
                data.append(str(getattr(self, field.lower())))
            
        out = '\t'.join(data)

        return out
    
    @classmethod
    def info_dict_to_str(cls, info):
        """Generate INFO string from dictionary
        
        Args:
            info: (dict) name=one of info_fields, value=must be str/int
        Returns:
            String of INFO data
        """
        data = []
        for field in cls.info_fields:
            if info.has_key(field):
                if cls.meta_info[field][1] == 'Flag':
                    data.append(field)
                else:
                    data.append('%s=%s' % (field, info[field]))
        return ';'.join(data)
    
    @classmethod
    def header(cls, source=None, reference_url=None, assembly_url=None):
        """Generate header lines"""
        lines = []
        line = '##fileformat=VCFv%s' % cls.version
        lines.append(line)
        
        # fileDate
        line = '##fileDate=%s' % time.strftime("%Y%m%d")
        lines.append(line)
        
        # source
        if source is not None:
            line = '##source=%s' % source
            lines.append(line)
                
        # reference url
        if reference_url is not None:
            line = '##reference=%s' % reference_url
            lines.append(line)
                
        # assembly url
        if assembly_url is not None:
            line = '##assembly=%s' % assembly_url
            lines.append(line)
            
        # INFO
        for item, info in cls.meta_info.items():
            line = '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">' % (item, info[0], info[1], info[2])
            lines.append(line)
            
        # ALT
        for item, desc in cls.meta_alt.items():
            line = '##ALT=<ID=%s,Description="%s">' % (item, desc)
            lines.append(line)
                        
        # columns
        line = '#%s' % '\t'.join(cls.fields)
        lines.append(line)
                    
        return '\n'.join(lines)

