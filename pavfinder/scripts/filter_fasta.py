import argparse
from collections import OrderedDict

def parse(fasta):
    seqs = OrderedDict()
    header = None
    seq = ''
    with open(fasta, 'r') as ff:
        for line in ff:
            if line[0] == '>':
                if header and seq:
                    seqs[header] = seq
                    header = None
                    seq = ''
                header = line.rstrip()
            else:
                seq += line.rstrip()
        if header and seq:
            seqs[header] = seq

    return seqs

def output(seqs, out_file, min_len):
    with open(out_file, 'w') as out:
        for header, seq in seqs.items():
            if len(seq) < min_len:
                continue
            out.write('{}\n{}\n'.format(header, seq))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=str, help="input fasta")
    parser.add_argument("min_len", type=int, help="minimum length")
    parser.add_argument("out", type=str, help="output fasta")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    seqs = parse(args.fasta)
    output(seqs, args.out, args.min_len)

if __name__ == '__main__':
    main()
