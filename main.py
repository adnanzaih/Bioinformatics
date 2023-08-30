from DNAToolkit import *

dna = "TAAAAGTGGGGCCAAAGTCTGGGCAGGCATGATGAT"

result = countNucFrequency(dna)
print(" ".join([str(val) for key, val in result.items()]))
#
# print(colored(dna))
# print(transcription(dna))
# print(reverse_complement(dna))
# print(gc_content(dna))
# print(gc_content_subset(dna, 5))
for _ in gen_reading_frames(dna):
    print(_)
