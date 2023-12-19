import itertools
from Bio import AlignIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner

'''
This represents a custom MSA identity and similarity script. I did not
create this by myself and was part of a collaboration with chatGPT.
The similarity scoring matrix is set to "BLOSUM62, but this can be changed.
The same is true for the gap scores. 
'''

def calculate_identity(seq1, seq2):
    matches = sum(s1 == s2 for s1, s2 in zip(seq1, seq2))
    return matches / len(seq1)

def calculate_similarity(seq1, seq2):
    blosum62 = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    valid_characters = set(blosum62.alphabet)
    filtered_seq1 = "".join(char for char in seq1 if char in valid_characters)
    filtered_seq2 = "".join(char for char in seq2 if char in valid_characters)

    alignments = aligner.align(filtered_seq1, filtered_seq2)
    alignment_score = alignments[0].score
    
    max_score_seq1 = sum(blosum62[char, char] for char in filtered_seq1)
    max_score_seq2 = sum(blosum62[char, char] for char in filtered_seq2)
    max_possible_score = min(max_score_seq1, max_score_seq2)

    return alignment_score / max_possible_score

#this is where you change the name and file path of your .clustal_num file
alignment = AlignIO.read("./outlier.aln", "clustal")

identities = []
similarities = []

for seq1, seq2 in itertools.combinations(alignment, 2):
    identity = calculate_identity(seq1.seq, seq2.seq)
    similarity = calculate_similarity(seq1.seq, seq2.seq)
    identities.append(identity)
    similarities.append(similarity)

average_identity = sum(identities) / len(identities)
average_similarity = sum(similarities) / len(similarities)

print("Average percent identity:", average_identity * 100)
print("Average percent similarity:", average_similarity * 100)


