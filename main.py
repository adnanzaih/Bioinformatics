from DNAToolkit import *

dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
result = countNucFrequency(dna)
print(" ".join([str(val) for key, val in result.items()]))