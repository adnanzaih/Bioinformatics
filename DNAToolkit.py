import collections
from structures import *
from utilities import *

def validateSeq(seq):
    """Check the sequence to make sure it is a valid DNA string"""
    tmpseq = seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

def countNucFrequency(seq):
    """Count nucleotides in a given sequence. Return a dictionary"""
    _return = collections.Counter(seq)
    return dict(sorted(_return.items(),key = lambda i: i[0]))

def transcription(seq):
    """DNA -> RNA Transcription. Replacing Thymine with Uracil"""
    return seq.replace("T", "U")

def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    #return ''.join([DNA_ReverseComplement[nuc] for nuc in dna_seq])[::-1]
    mapping = str.maketrans("ATCG", 'TAGC')
    return seq.translate(mapping)[::-1]

def gc_content(seq):
    """GC Content in a DNA/RNA sequence"""
    return round((seq.count("C") + seq.count("G")) / len(seq) * 100, 6)

def gc_content_subset(seq, k=20):
    """GC Content in a DNA/RNA sub-sequence length k, k = 20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res

def translate_seq(seq, init_pos=0):
    "Translate a DNA sequence into an aminoacid sequence"
    return [DNA_Codons[seq[pos:pos+3]] for pos in range(init_pos, len(seq)-2,3)]

def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacide in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq)-2, 3):
        if DNA_Codons[seq[i:i+3]] == aminoacid:
            tmpList.append(seq[i:i+3])

    freqDict = dict(collections.Counter(tmpList))
    totalWeight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] /  totalWeight, 2)
    return freqDict

def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including the reverse complement"""
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames