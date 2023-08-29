import collections

Nucleotides = {
    "A",
    "C",
    "G",
    "T"
}

def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

def countNucFrequency(dna_seq):
    _return = collections.Counter(dna_seq)
    return dict(sorted(_return.items(),key = lambda i: i[0]))